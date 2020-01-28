#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Copyright 2018-2019 (C) DOE JGI, LBL
# Seung-Jin Sul (ssul@lbl.gov)

import os
import sys
import unittest
import socket

# For unittest
ROOT_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir)
LIB_DIR = '%s' % ROOT_DIR
if LIB_DIR not in sys.path:
    sys.path.append(LIB_DIR)

from .RabbitmqConnection import *

DEFAULT_POOL = ["small", "medium", "large", "xlarge"]

class JtmInterface(object):
    """
    Class for jtm-* CLI tools
    """
    def __init__(self, taskType):
        self.taskType = taskType
        self.rabbitConnection = RmqConnectionHB()
        self.connection = self.rabbitConnection.open()
        self.channel = self.connection.channel()
        self.channel.exchange_declare(exchange=JGI_JTM_MAIN_EXCH,
                                      exchange_type='direct',
                                      durable=True,
                                      auto_delete=False)

        self.response = None

        # OLD
        # result = self.channel.queue_declare('_jtm_submit.' + CNAME,

        # NEW
        # Create a temp random queue
        result = self.channel.queue_declare('',
                                            durable=True,
                                            # exclusive=True,
                                            exclusive=False,
                                            auto_delete=True)
        self.callback_queue = result.method.queue
        self.channel.queue_bind(exchange=JGI_JTM_MAIN_EXCH,
                                queue=self.callback_queue,
                                routing_key=self.callback_queue)
        self.channel.basic_qos(prefetch_count=1)

        if int(PIKA_VER[0]) < 1:  # v0.13.1
            self.channel.basic_consume(self.on_response,
                                       queue=self.callback_queue,
                                       no_ack=False)
        else:  # v1.0.1 or higher
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
        jsonStr = ""
        jtmHostName = ""

        if "task_file" in kw and kw["task_file"]:
            with open(kw["task_file"], "r") as jf:
                for l in jf:
                    if l:
                        jsonStr += l
        elif "task_json" in kw and kw["task_json"]:
            jsonStr = '{"command": "%s"}' % kw["task_json"]

        if "task_id" in kw and kw["task_id"] == 0:  # jtm-submit
            try:
                jsonData = json.loads(jsonStr)
            except ValueError:
                print(jsonStr)
                raise
        else:  # jtm-status, jtm-kill, jtm-resource, jtm-check-*
            jsonData = {}

        # pprint.pprint(jsonData)

        if "task_pool" in kw and kw["task_pool"]:
            jsonData["task_pool"] = kw["task_pool"]
        if "jtm_host_name" in kw and kw["jtm_host_name"]:
            jtmHostName = jsonData["jtm_host_name"] = kw["jtm_host_name"].replace(".", "_")

        jsonData['task_type'] = self.taskType
        jsonData['task_id'] = kw["task_id"] if "task_id" in kw else 0

        # Check if exlusive option is set and
        # Pass Cromwell job id and step name to jtm manager
        excluTaskPoolNamePostfix = ""
        if ("shared" in kw and kw["shared"] == 0) or \
           ("pool" in jsonData and "shared" in jsonData["pool"] and jsonData["pool"]["shared"] == 0):
            # Note: here "job_id" is a Cromwell workflow job name sent by "jtm-submit -jid ${job_name}"
            #  from cromwell.conf
            # The "shared" means the dynamic pool won't be shared among multiple workflows
            # and it will be disposed after this wf is completed.
            # If shared=0, append cromwell job id to the user pool name
            # so that the pool can be removed once the wf is done.
            # ex) cromwell_fe010880_stepA ==> extract fe010880
            # Ref) https://issues.jgi-psf.org/browse/JAWS-8
            #
            if "pool" in jsonData and "job_id" in jsonData["pool"] and jsonData["pool"]["job_id"]:
                excluTaskPoolNamePostfix = jsonData["pool"]["job_id"].split('_')[1]
                jsonData["pool"]["name"] = jsonData["pool"]["name"] + '_' + excluTaskPoolNamePostfix
            elif "job_id" in kw and kw["job_id"]:
                jsonData["job_id"] = kw["job_id"]
                jsonData["shared"] = int(kw["shared"])
                excluTaskPoolNamePostfix = kw["job_id"].split('_')[1]
                kw['pool_name'] = kw['pool_name'] + '_' + excluTaskPoolNamePostfix

        # If cromwell task and custom pool setting found, create a separate pool for the tasks!
        if "job_time" in kw:
            if kw["job_time"] and kw["node_mem"] and kw["num_core"] and \
               kw["pool_name"] and kw["job_time"] != "" and kw["node_mem"] != "" and \
               kw["num_core"] != 0 and kw["pool_name"] != "" and kw["node"] != 0 and \
               kw["job_time"] != "00:00:00" and kw["node_mem"] != "0G" and kw["num_core"] != 0 and \
               kw["pool_name"] != "default" and kw["pool_name"] not in DEFAULT_POOL:
                jsonData["pool"] = {}
                jsonData["pool"]["time"] = kw["job_time"]
                jsonData["pool"]["cpu"] = kw["num_core"]
                jsonData["pool"]["mem"] = kw["node_mem"]
                jsonData["pool"]["name"] = kw["pool_name"]
                jsonData["pool"]["cluster"] = kw["jtm_host_name"]
                jsonData["pool"]["nwpn"] = kw["nwpn"] if 'nwpn' in kw else 1  # number of workers per node
                jsonData["pool"]["node"] = kw["node"] if 'node' in kw else 1  # number of nodes
                jsonData["pool"]["shared"] = int(kw["shared"])
                jsonData["pool"]["constraint"] = kw["constraint"]
                jsonData["pool"]["qos"] = kw["qos"]
                jsonData["pool"]["account"] = kw["account"]

        # For the command like, "jtm-submit -cr 'ls' -cl cori -p test"
        if "pool" in jsonData and "name" in jsonData["pool"] and jsonData["pool"]["name"] is not None:
            pass
        else:
            if "pool_name" in kw and kw["pool_name"] and kw["pool_name"] != "default":
                jsonData["pool"] = {}
                jsonData["pool"]["name"] = kw['pool_name']


        msgZipped = zdumps(json.dumps(jsonData))


        # Determine the destination queues

        # USER_NAME = getpass.getuser()
        userName = USER_NAME  # Config.py. For now, it's fixed as "jaws"

        # If jtmHostName is not set, try to find it
        if not jtmHostName:
            if "pool" in jsonData and jsonData["pool"] and "cluster" in jsonData["pool"]:
                jtmHostName = jsonData["pool"]["cluster"]
                # if "name" in jsonData["pool"]:
                #     customPoolName = jsonData["pool"]["name"]
            elif "jtm_host_name" in jsonData and jsonData["jtm_host_name"]:
                jtmHostName = jsonData["jtm_host_name"]
            else:
                if "NERSC_HOST" in os.environ:
                    jtmHostName = os.environ["NERSC_HOST"]
                elif "JTM_HOST_NAME" in os.environ:  # for custom name like ["aws' | 'olcf' | 'pc']
                    jtmHostName = os.environ["JTM_HOST_NAME"]
                elif "HOSTNAME" in os.environ:
                    jtmHostName = os.environ["HOSTNAME"]
                else:
                    jtmHostName = socket.gethostname()

        # Prepare rmq message
        jtmHostName = jtmHostName.replace(".", "_")
        cName = jtmHostName + "." + userName  # username = jtm

        # It's not good to have here again but it's for dealing with multiple jtm instances
        jtmTaskRequestQ = "_jtm_task_request_queue" + "." + cName
        # jtmTaskResultQ = "_jtm_task_result_queue" + "." + cName + "." + self.taskType

        # Note: declare and bind are needed to keep jtm-submit requests in the queue
        #  and to let the manager consume it
        self.channel.queue_declare(queue=JTM_TASK_REQUEST_Q,
                                   durable=True,
                                   exclusive=False,
                                   auto_delete=True)
        self.channel.queue_bind(exchange=JGI_JTM_MAIN_EXCH,
                                queue=JTM_TASK_REQUEST_Q,
                                routing_key=JTM_TASK_REQUEST_Q)

        if "log_level" in kw and kw["log_level"] == "debug":
            print("kw")
            pprint.pprint(kw)
            print("jsonData")
            pprint.pprint(jsonData)
            print("jtmTaskRequestQ = %s" % jtmTaskRequestQ)
            # print "jtmTaskResultQ =", jtmTaskResultQ

        self.corr_id = str(uuid.uuid4())
        # self.channel.confirm_delivery()  # 11192018 to test task is not discarded
        self.channel.basic_publish(exchange=JGI_JTM_MAIN_EXCH,
                                   routing_key=jtmTaskRequestQ,
                                   properties=pika.BasicProperties(
                                       delivery_mode=2,  # added for fixing jtm-submit timeout
                                       reply_to=self.callback_queue,
                                       correlation_id=self.corr_id),
                                   body=msgZipped)

        cnt = 0
        try:
            while self.response is None:
                # time_limit is actually acting as interval in secs
                # Note: suggested upper bound on processing time in seconds. The actual blocking time depends on the
                #  granularity of the underlying ioloop. Zero means return as soon as possible.
                #  None means there is no limit on processing time
                #  and the functionwill block until I/O produces actionable events.
                #
                # Fixme: timeout. time_limit=0 or time_limit=None doesn't work. Set long enough processing time
                #  fixed it! 04162019
                #
                # time.sleep(1)
                # self.connection.process_data_events(time_limit=JTMINTERFACE_TIMEOUT)
                self.connection.process_data_events(time_limit=1.0)

                cnt += 1
                if cnt == JTMINTERFACE_MAX_TRIAL:
                    make_dir_p(os.path.join(JTM_LOG, "jtm-submit"))
                    eprint("Failed to get a reply from the manager: {} {}".format(jsonData, self.response))
                    with open(os.path.join(JTM_LOG, "jtm-submit", "jtm_submit_%s" % (datetime.datetime.now().strftime("%Y-%m-%d"))), 'a') as jslogf:
                        jslogf.write("{} {}\n".format(jsonData, self.response))

                    if self.response is None:  # still none?
                        self.response = -88
                    break
        except Exception as e:
            print("no task id returned")
            print(e)


        return self.response


    def close(self):
        self.channel.close()
        self.connection.close()


class TestJtmInterface(unittest.TestCase):
    def testTask(self):
        print(" [x] Submitting a task")
        jtm_submit = JtmInterface('task')
        response = jtm_submit.call(task_file="../example/jtm_task_simple.json")
        self.assertIsNotNone(response)
        print(" [.] Got %r\n" % response)

    # def testStatus(self):
    #     print(" [x] Requesting status")
    #     jtm_status = JtmInterface('status')
    #     response = jtm_status.call(task_id=1)
    #     self.assertIsNotNone(response)
    #     print(" [.] Got %r\n" % response)
    #
    # def testKill(self):
    #     print(" [x] Cancelling a task")
    #     jtm_kill = JtmInterface('')
    #     response = jtm_kill.call(task_id=1)
    #     self.assertIsNotNone(response)
    #     print(" [.] Got %r\n" % response)


# Use this for testing purpose
if __name__ == "__main__":
    unittest.main()
