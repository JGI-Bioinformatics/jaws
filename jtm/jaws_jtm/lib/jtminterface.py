#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)
import os
import pprint
import shortuuid
import json
import uuid
import datetime
import pika

from jaws_jtm.lib.rabbitmqconnection import RmqConnectionHB
from jaws_jtm.lib.msgcompress import zdumps, zloads
from jaws_jtm.lib.run import make_dir, eprint


class JtmInterface(object):
    """
    Class for jtm-* CLI tools
    """
    def __init__(self, task_type, ctx=None, info_tag=None,):
        self.task_type = task_type
        self.config = ctx.obj['config']
        self.rmq_conn = RmqConnectionHB(config=self.config)
        self.connection = self.rmq_conn.open()
        self.channel = self.connection.channel()
        self.response = None
        self.JTM_LOG = self.config.configparser.get("JTM", "log_dir")
        self.JGI_JTM_MAIN_EXCH = self.config.configparser.get("JTM", "jgi_jtm_main_exch")
        self.USER_NAME = self.config.configparser.get("SITE", "user_name")
        self.JTM_TASK_REQUEST_Q = self.config.configparser.get("JTM", "jtm_task_request_q")
        self.JTMINTERFACE_MAX_TRIAL = self.config.configparser.getint("JTM", "jtminterface_max_trial")
        self.JTM_HOST_NAME = self.config.configparser.get("SITE", "jtm_host_name")
        self.DEBUG = ctx.obj['debug']
        self.corr_id = str(uuid.uuid4())

        self.channel.exchange_declare(exchange=self.JGI_JTM_MAIN_EXCH,
                                      exchange_type='direct',
                                      durable=True,
                                      auto_delete=False)

        # OLD
        # result = self.channel.queue_declare('_jtm_submit.' + CNAME,

        # NEW
        # Create a temp random queue
        if info_tag:
            uniq_queue_name = "_jtm_client_" + '_'.join([str(shortuuid.uuid()), task_type, str(info_tag)])
            self.channel.queue_declare(uniq_queue_name,
                                       durable=True,
                                       # exclusive=True,
                                       exclusive=False,
                                       auto_delete=True)
            self.callback_queue = uniq_queue_name
        else:
            result = self.channel.queue_declare('',
                                                durable=True,
                                                # exclusive=True,
                                                exclusive=False,
                                                auto_delete=True)
            self.callback_queue = result.method.queue

        # print(self.callback_queue)

        self.channel.queue_bind(exchange=self.JGI_JTM_MAIN_EXCH,
                                queue=self.callback_queue,
                                routing_key=self.callback_queue)
        self.channel.basic_qos(prefetch_count=1)

        self.channel.basic_consume(queue=self.callback_queue,
                                   on_message_callback=self.on_response,
                                   auto_ack=False)

    def on_response(self, ch, method, props, body):
        if self.corr_id == props.correlation_id:
            self.response = zloads(body) if body else None
            self.channel.basic_ack(delivery_tag=method.delivery_tag)
        else:
            # self.channel.basic_nack(delivery_tag=method.delivery_tag)
            # If corr_id is not for me, reject and requeue it
            self.channel.basic_reject(delivery_tag=method.delivery_tag, requeue=True)

    def call(self, **kw):
        """
        Handles functions for CLI
        :param kw:
        :return:
        """
        json_str = ""
        jtm_host_name = ""

        if "task_file" in kw and kw["task_file"]:
            with open(kw["task_file"], "r") as jf:
                for task_file_read_line in jf:
                    if task_file_read_line:
                        json_str += task_file_read_line
        elif "task_json" in kw and kw["task_json"]:
            json_str = '{"command": "%s"}' % kw["task_json"]

        if "task_id" in kw and kw["task_id"] == 0:  # jtm-submit
            try:
                json_data_dict = json.loads(json_str)
            except ValueError:
                print(json_str)
                raise
        else:  # jtm-status, jtm-kill, jtm-resource, jtm-check-*
            json_data_dict = {}

        if "task_pool" in kw and kw["task_pool"]:
            json_data_dict["task_pool"] = kw["task_pool"]
        if "jtm_host_name" in kw and kw["jtm_host_name"]:
            jtm_host_name = json_data_dict["jtm_host_name"] = kw["jtm_host_name"].replace(".", "_")

        json_data_dict['task_type'] = self.task_type
        json_data_dict['task_id'] = kw["task_id"] if "task_id" in kw else 0
        json_data_dict['slurm_info'] = True if "slurm_info" in kw else False

        # Check if exlusive option is set (shared=0) and
        # Pass Cromwell job id and step name to jtm manager
        exclusive_task_pool_name_postfix = ""
        if ("shared" in kw and kw["shared"] == 0) or \
           ("pool" in json_data_dict and "shared" in json_data_dict["pool"] and json_data_dict["pool"]["shared"] == 0):
            # Note: here "job_id" is a Cromwell workflow job name sent by "jtm-submit -jid ${job_name}"
            #  from cromwell.conf
            # The "shared" means the dynamic pool won't be shared among multiple workflows
            # and it will be disposed after this wf is completed.
            # If shared=0, append cromwell job id to the user pool name
            # so that the pool can be removed once the wf is done.
            # ex) cromwell_fe010880_stepA ==> extract fe010880
            # Ref) https://issues.jgi-psf.org/browse/JAWS-8
            #
            if "pool" in json_data_dict and "job_id" in json_data_dict["pool"] and json_data_dict["pool"]["job_id"]:
                exclusive_task_pool_name_postfix = json_data_dict["pool"]["job_id"].split('_')[1]
                json_data_dict["pool"]["name"] = json_data_dict["pool"]["name"] + '_' + exclusive_task_pool_name_postfix
            elif "job_id" in kw and kw["job_id"]:
                json_data_dict["job_id"] = kw["job_id"]
                json_data_dict["shared"] = int(kw["shared"])
                exclusive_task_pool_name_postfix = kw["job_id"].split('_')[1]
                kw['pool_name'] = kw['pool_name'] + '_' + exclusive_task_pool_name_postfix

        # If cromwell task and custom pool setting found, create a separate pool for the tasks!
        if "job_time" in kw:
            if kw["job_time"] and kw["node_mem"] and kw["num_core"] and \
               kw["pool_name"] and kw["job_time"] != "" and kw["node_mem"] != "" and \
               kw["num_core"] != 0 and kw["pool_name"] != "" and kw["node"] != 0 and \
               kw["job_time"] != "00:00:00" and kw["node_mem"] != "0G" and kw["num_core"] != 0:
                json_data_dict["pool"] = {}
                json_data_dict["pool"]["time"] = kw["job_time"]
                json_data_dict["pool"]["cpu"] = kw["num_core"]
                json_data_dict["pool"]["mem"] = kw["node_mem"]
                json_data_dict["pool"]["name"] = kw["pool_name"]
                json_data_dict["pool"]["cluster"] = kw["jtm_host_name"]
                json_data_dict["pool"]["nwpn"] = kw["nwpn"] if 'nwpn' in kw else 1  # number of workers per node
                json_data_dict["pool"]["node"] = kw["node"] if 'node' in kw else 1  # number of nodes
                json_data_dict["pool"]["shared"] = int(kw["shared"])
                json_data_dict["pool"]["constraint"] = kw["constraint"]
                json_data_dict["pool"]["qos"] = kw["qos"]
                json_data_dict["pool"]["account"] = kw["account"]

        # For the command like, "jtm-submit -cr 'ls' -cl cori -p test"
        if "pool" in json_data_dict and "name" in json_data_dict["pool"] and json_data_dict["pool"]["name"] is not None:
            pass
        else:
            if "pool_name" in kw and kw["pool_name"] and kw["pool_name"] != "default":
                json_data_dict["pool"] = {}
                json_data_dict["pool"]["name"] = kw['pool_name']

        msg_zipped = zdumps(json.dumps(json_data_dict))

        # Determine the destination queues

        # self.USER_NAME = getpass.getuser()
        user_name = self.USER_NAME  # Config.py. For now, it's fixed as "jaws"
        jtm_host_name = self.JTM_HOST_NAME

        # If jtm_host_name is not set, try to find it
        # if not jtm_host_name:
        #     if "pool" in json_data_dict and json_data_dict["pool"] and "cluster" in json_data_dict["pool"]:
        #         jtm_host_name = json_data_dict["pool"]["cluster"]
        #         # if "name" in json_data_dict["pool"]:
        #         #     customPoolName = json_data_dict["pool"]["name"]
        #     elif "jtm_host_name" in json_data_dict and json_data_dict["jtm_host_name"]:
        #         jtm_host_name = json_data_dict["jtm_host_name"]
        #     else:
        #         if "NERSC_HOST" in os.environ:
        #             jtm_host_name = os.environ["NERSC_HOST"]
        #         elif "self.JTM_HOST_NAME" in os.environ:  # for custom name like ["aws' | 'olcf' | 'pc']
        #             jtm_host_name = os.environ["self.JTM_HOST_NAME"]
        #         elif "HOSTNAME" in os.environ:
        #             jtm_host_name = os.environ["HOSTNAME"]
        #         else:
        #             jtm_host_name = socket.gethostname()

        # Prepare rmq message
        jtm_host_name = jtm_host_name.replace(".", "_")
        host_and_user_name = jtm_host_name + "." + user_name

        # It's not good to have here again but it's for dealing with multiple jtm instances
        jtmTaskRequestQ = "_jtm_task_request_queue" + "." + host_and_user_name
        # jtmTaskResultQ = "_jtm_task_result_queue" + "." + host_and_user_name + "." + self.task_type

        # Note: declare and bind are needed to keep jtm-submit requests in the queue
        #  and to let the manager consume it
        self.channel.queue_declare(queue=self.JTM_TASK_REQUEST_Q,
                                   durable=True,
                                   exclusive=False,
                                   auto_delete=True)
        self.channel.queue_bind(exchange=self.JGI_JTM_MAIN_EXCH,
                                queue=self.JTM_TASK_REQUEST_Q,
                                routing_key=self.JTM_TASK_REQUEST_Q)

        if self.DEBUG:
            print("kw")
            pprint.pprint(kw)
            print("json_data_dict")
            pprint.pprint(json_data_dict)
            print(("jtmTaskRequestQ = %s" % jtmTaskRequestQ))
            # print "jtmTaskResultQ =", jtmTaskResultQ

        # self.channel.confirm_delivery()  # 11192018 to test task is not discarded
        self.channel.basic_publish(exchange=self.JGI_JTM_MAIN_EXCH,
                                   routing_key=jtmTaskRequestQ,
                                   properties=pika.BasicProperties(
                                       delivery_mode=2,  # added for fixing jtm-submit timeout
                                       reply_to=self.callback_queue,
                                       correlation_id=self.corr_id),
                                   body=msg_zipped)

        cnt = 0
        try:
            while self.response is None:
                # time_limit is actually acting as interval in secs
                # Note: suggested upper bound on processing time in seconds. The actual blocking time depends on the
                #  granularity of the underlying ioloop. Zero means return as soon as possible.
                #  None means there is no limit on processing time
                #  and the function will block until I/O produces actionable events.
                #
                # Fixme: timeout. time_limit=0 or time_limit=None doesn't work. Set long enough processing time
                #  fixed it! 04162019
                #
                # time.sleep(1)
                # self.connection.process_data_events(time_limit=JTMINTERFACE_TIMEOUT)
                self.connection.process_data_events(time_limit=1.0)

                cnt += 1
                if cnt == self.JTMINTERFACE_MAX_TRIAL:
                    make_dir(os.path.join(self.JTM_LOG, "jtm-submit"))
                    eprint("Failed to get a reply from the manager: {} {}".format(json_data_dict, self.response))
                    jtm_submit_log_file = os.path.join(self.JTM_LOG, "jtm-submit", "jtm_submit_%s"
                                                       % (datetime.datetime.now().strftime("%Y-%m-%d")))
                    with open(jtm_submit_log_file, 'a') as jslogf:
                        jslogf.write("{} {}\n".format(json_data_dict, self.response))
                    os.chmod(jtm_submit_log_file, 0o777)
                    if self.response is None:  # still none?
                        self.response = -88
                    break
        except Exception as e:
            eprint("No task id returned")
            eprint(e)

        return self.response

    def close(self):
        self.channel.close()
        self.connection.close()
