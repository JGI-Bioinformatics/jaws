#! /usr/bin/env python
# pylint: disable=C0111,C0103,R0205
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)

"""

JGI task managemr

Example of task processing scenario
1. jtm-submiit sends a msg to "jgi_main_exchange" with "jtm_task_request_queue" tag.
2. jtm-manager listens to "jtm_task_request_queue" which is bound to "jgi_main_exchange"
3. When a task is published, jtm-manager takes it and sends it to a pool
   (to jgi_jtm_inner_main_exchange)
4. Workers listen to jgi_jtm_inner_request_queue which is bound to "jgi_jtm_inner_main_exchange"
5. When a task is completed, a worker sends a result msg to "jgi_jtm_inner_main_exchange" with
   "jgi_jtm_inner_result_queue" tag
6. jtm-manager listens to "jgi_jtm_inner_result_queue" queue. When a result is ready,
   takes and updates tables

"""
import sys
import multiprocessing as mp
import time
import json
from math import ceil
import uuid
import datetime
import shortuuid
import os
import signal
import re
import amqpstorm
import psutil

from jaws_jtm.common import setup_custom_logger, logger
from jaws_jtm.lib.sqlstmt import JTM_SQL
from jaws_jtm.lib.rabbitmqconnection import RmqConnectionAmqpstorm, JtmAmqpstormBase
from jaws_jtm.lib.dbutils import DbSqlMysql
from jaws_jtm.lib.run import (
    pad_string_path,
    make_dir,
    run_sh_command,
    extract_cromwell_id,
)
from jaws_jtm.lib.msgcompress import zdumps, zloads
from jaws_rpc import rpc_client, rpc_server, responses


# --------------------------------------------------------------------------------------------------
# Globals
# --------------------------------------------------------------------------------------------------
NUM_TOTAL_WORKERS = mp.Value("i", 0)


class WorkerResultReceiver(JtmAmqpstormBase):
    """
    Process result from worker

    """

    def start(self):
        """Start the WorkerResultReceiver.
        :return:
        """
        if not self.connection:
            self.create_connection()
        while True:
            try:
                channel = self.connection.channel()
                channel.exchange.declare(
                    exchange=self.jtm_inner_main_exch,
                    exchange_type="direct",
                    durable=True,
                    auto_delete=False,
                )
                channel.queue.declare(
                    queue=self.inner_result_queue_name,
                    durable=True,
                    exclusive=False,
                    auto_delete=True,
                )
                channel.queue.bind(
                    exchange=self.jtm_inner_main_exch,
                    queue=self.inner_result_queue_name,
                    routing_key=self.inner_result_queue_name,
                )
                channel.basic.consume(
                    self.process_result, self.inner_result_queue_name, no_ack=False
                )
                channel.basic.qos(prefetch_count=1)
                channel.start_consuming()
                if not channel.consumer_tags:
                    channel.close()
            except amqpstorm.AMQPError as why:
                logger.exception(why)
                self.create_connection()
            except KeyboardInterrupt:
                self.connection.close()
                break

    def process_result(self, message):
        """
        Callback for handling request from JtmInterface

        :param message:
        :return:
        """
        msg_unzipped = json.loads(zloads(message.body))
        logger.info("New result: {}".format(msg_unzipped))
        assert "task_id" in msg_unzipped, "task ID not found from the result package"
        task_id = int(msg_unzipped["task_id"])
        assert task_id >= 0, "Found invalid task ID"
        done_flag = (
            int(msg_unzipped["done_flag"]) if "done_flag" in msg_unzipped else -2
        )
        ret_msg = msg_unzipped.get("ret_msg", "")
        a_worker_id = msg_unzipped.get("worker_id", "")
        host_name = msg_unzipped.get("host_name", "")

        if ret_msg != "hb":
            db = DbSqlMysql(config=self.config)
            db.execute(
                JTM_SQL["update_tasks_doneflag_by_taskid"]
                % dict(task_id=task_id, done_flag=done_flag),
                debug=DEBUG,
            )
            db.commit()
            db.close()

            # Check the task status code
            if done_flag > 0:
                task_status_int = self.task_status["success"]  # 4
            else:
                if done_flag == -1:
                    task_status_int = self.task_status["outputerror"]
                elif done_flag == -2:
                    task_status_int = self.task_status["failed"]
                elif done_flag == -3:
                    task_status_int = self.task_status["outofresource"]
                elif done_flag == -4:
                    task_status_int = self.task_status["terminated"]
                elif done_flag == -5:
                    task_status_int = self.task_status["failed"]
                elif done_flag == -6:
                    task_status_int = self.task_status["timeout"]
                elif done_flag == -7:
                    task_status_int = self.task_status["lostconnection"]
                else:
                    logger.critical("Unknown return code {}".format(done_flag))
                    raise OSError(2)

            db = DbSqlMysql(config=self.config)
            status_now = db.selectScalar(
                JTM_SQL["select_status_runs_by_taskid"] % dict(task_id=task_id),
                debug=False,
            )
            logger.debug(f"status now {task_id} = {status_now}")

            if not STANDALONE:
                logger.debug(
                    f"status change msg {task_id}: {status_now} => succeed/failed"
                )
                if status_now in (
                    self.task_status["ready"],
                    self.task_status["queued"],
                    self.task_status["pending"],
                    self.task_status["running"],
                ):
                    ret_status = 4 if task_status_int > 0 else -2
                    if status_now == self.task_status["ready"]:
                        logger.debug("process result --> ready to queued")
                        send_update_task_status_msg(
                            task_id, TASK_STATUS["ready"], TASK_STATUS["queued"]
                        )
                        db.execute(
                            JTM_SQL["update_runs_status_by_taskid"]
                            % dict(
                                status_id=TASK_STATUS["queued"],
                                task_id=task_id,
                                now=time.strftime("%Y-%m-%d %H:%M:%S"),
                            ),
                            debug=DEBUG,
                        )
                        db.commit()
                        logger.debug("process result --> queued to pending")
                        send_update_task_status_msg(
                            task_id, TASK_STATUS["queued"], TASK_STATUS["pending"]
                        )
                        db.execute(
                            JTM_SQL["update_runs_status_by_taskid"]
                            % dict(
                                status_id=TASK_STATUS["pending"],
                                task_id=task_id,
                                now=time.strftime("%Y-%m-%d %H:%M:%S"),
                            ),
                            debug=DEBUG,
                        )
                        db.commit()
                        logger.debug("process result --> pending to running")
                        send_update_task_status_msg(
                            task_id, TASK_STATUS["pending"], TASK_STATUS["running"]
                        )
                        db.execute(
                            JTM_SQL["update_runs_status_by_taskid"]
                            % dict(
                                status_id=TASK_STATUS["running"],
                                task_id=task_id,
                                now=time.strftime("%Y-%m-%d %H:%M:%S"),
                            ),
                            debug=DEBUG,
                        )
                        db.commit()
                    elif status_now == self.task_status["queued"]:
                        logger.debug("process result --> queued to pending")
                        send_update_task_status_msg(
                            task_id, TASK_STATUS["queued"], TASK_STATUS["pending"]
                        )
                        db.execute(
                            JTM_SQL["update_runs_status_by_taskid"]
                            % dict(
                                status_id=TASK_STATUS["pending"],
                                task_id=task_id,
                                now=time.strftime("%Y-%m-%d %H:%M:%S"),
                            ),
                            debug=DEBUG,
                        )
                        db.commit()
                        logger.debug("process result --> pending to running")
                        send_update_task_status_msg(
                            task_id, TASK_STATUS["pending"], TASK_STATUS["running"]
                        )
                        db.execute(
                            JTM_SQL["update_runs_status_by_taskid"]
                            % dict(
                                status_id=TASK_STATUS["running"],
                                task_id=task_id,
                                now=time.strftime("%Y-%m-%d %H:%M:%S"),
                            ),
                            debug=DEBUG,
                        )
                        db.commit()
                    elif status_now == self.task_status["pending"]:
                        logger.debug("process result --> pending to running")
                        send_update_task_status_msg(
                            task_id, TASK_STATUS["pending"], TASK_STATUS["running"]
                        )
                        db.execute(
                            JTM_SQL["update_runs_status_by_taskid"]
                            % dict(
                                status_id=TASK_STATUS["running"],
                                task_id=task_id,
                                now=time.strftime("%Y-%m-%d %H:%M:%S"),
                            ),
                            debug=DEBUG,
                        )
                        db.commit()

                    logger.debug("process result --> running to succeeded/failed")
                    send_update_task_status_msg(
                        task_id,
                        TASK_STATUS["running"],
                        ret_status,
                        fail_code=done_flag if done_flag < 0 else None,
                        reason=ret_msg,
                    )

            db.execute(
                JTM_SQL["update_tasks_doneflag_by_taskid"]
                % dict(task_id=task_id, done_flag=done_flag),
                debug=DEBUG,
            )
            db.commit()
            db.close()

            time.sleep(CONFIG.configparser.getfloat("JTM", "result_receive_interval"))

            # Sometimes workerId2 ==> 0
            # so wait a little bit if that happened
            a_worker_id_to_check = 0
            db = DbSqlMysql(config=CONFIG)
            rows = db.selectAll(
                JTM_SQL["select_workerid2_workers_by_wid"] % dict(worker_id=a_worker_id)
            )
            db.close()

            try:
                a_worker_id_to_check = int(rows[0][0])
            except Exception:
                a_worker_id_to_check = 0
            logger.debug(
                f"WorkerId2 from workerid = {a_worker_id_to_check}, {a_worker_id}"
            )

            # Just in case
            # If a worker has already been terminated, this worker might be deleted already.
            # Let's try only (wait_count * "worker_info_update_wait") seconds.
            wait_count = 0
            while a_worker_id_to_check == 0:
                wait_count += 1
                db = DbSqlMysql(config=CONFIG)
                rows = db.selectAll(
                    JTM_SQL["select_workerid2_workers_by_wid"]
                    % dict(worker_id=a_worker_id)
                )
                db.close()

                try:
                    a_worker_id_to_check = int(rows[0][0])
                except Exception:
                    a_worker_id_to_check = 0

                logger.debug(f"Try to get the real worker ID for {a_worker_id}")

                # Wait 20 times with 0.5 sec sleep
                time.sleep(
                    CONFIG.configparser.getfloat("JTM", "worker_info_update_wait")
                )
                if wait_count == 20:
                    a_worker_id_to_check = 0
                    break

            success = False
            while success is not True:
                success = True
                try:
                    db = DbSqlMysql(config=CONFIG)
                    db.execute(
                        JTM_SQL["update_runs_status_workerid2_by_taskid_2"]
                        % dict(
                            status_id=task_status_int,
                            wid2=a_worker_id_to_check,
                            now=time.strftime("%Y-%m-%d %H:%M:%S"),
                            task_id=task_id,
                        ),
                        debug=DEBUG,
                    )
                    db.commit()
                    db.close()
                    logger.debug(
                        "status {}, worker id {}, taskid {}".format(
                            task_status_int, a_worker_id_to_check, task_id
                        )
                    )
                except Exception as e:
                    logger.critical(e)
                    logger.critical(
                        "Failed to update runs table for status and workerid2."
                    )
                    logger.debug("Retry to update runs table for status and workerid2.")
                    success = False
                    time.sleep(
                        CONFIG.configparser.getfloat("JTM", "runs_info_update_wait")
                    )

            # Print report
            if done_flag == self.done_flag["success"]:  # 1
                logger.info(
                    "Task %s --> Success on worker/host, %s/%s",
                    task_id,
                    a_worker_id,
                    host_name,
                )
            elif (
                done_flag == self.done_flag["success with correct output file(s)"]
            ):  # 2
                logger.info(
                    "Task %s --> Success with valid output(s) on worker/host, %s/%s",
                    task_id,
                    a_worker_id,
                    host_name,
                )
            elif done_flag == self.done_flag["failed to check output file(s)"]:  # -1
                logger.info(
                    "Task %s --> %s, worker/host: %s/%s",
                    task_id,
                    ret_msg,
                    a_worker_id,
                    host_name,
                )
            elif done_flag == self.done_flag["failed to run user command"]:  # -2
                logger.info(
                    "Task %s --> Failed with non-zero exit code. stdout = %s, worker/host: %s/%s",
                    task_id,
                    ret_msg,
                    a_worker_id,
                    host_name,
                )
            elif done_flag == self.done_flag["failed with out-of-mem"]:  # -3
                logger.info(
                    "Task %s --> Failed with out-of-mem. stdout = %s, worker/host: %s/%s",
                    task_id,
                    ret_msg,
                    a_worker_id,
                    host_name,
                )
            elif done_flag == self.done_flag["failed with user termination"]:  # -4
                logger.info(
                    "Task %s --> Failed by user termination. stdout = %s, worker/host: %s/%s",
                    task_id,
                    ret_msg,
                    a_worker_id,
                    host_name,
                )
            elif (
                done_flag
                == self.done_flag["failed with input file or command not found"]
            ):  # -5
                logger.info(
                    "Task %s --> Failed with input file or command not found. stdout = %s, worker/host: %s/%s",
                    task_id,
                    ret_msg,
                    a_worker_id,
                    host_name,
                )
            elif done_flag == self.done_flag["failed with timeout"]:  # -6
                logger.info(
                    "Task %s --> Failed with timeout. stdout = %s, worker/host: %s/%s",
                    task_id,
                    ret_msg,
                    a_worker_id,
                    host_name,
                )
            else:
                logger.warning("Cannot recognize the return code: %d" % done_flag)
                logger.info(
                    "Task %s --> worker/host: %s/%s", task_id, a_worker_id, host_name
                )
                raise OSError(2, "Cannot recognize the DONE_FLAG return code")
        else:
            # If it is not a result msg, return it back to the exchange
            message.nack()

        # After this the result message will be deleted from RabbitMQ
        # If this worker crashes while running a user command, this task will
        # be sent to other workers available
        message.ack()


class JtmCommandRunner(JtmAmqpstormBase):
    """
    Process request from JtmInterface <- JTM Client

    """

    def start(self):
        """Start the JtmCommandRunners.
        :return:
        """
        if not self.connection:
            self.create_connection()
        while True:
            try:
                channel = self.connection.channel()
                channel.exchange.declare(
                    exchange=self.jgi_jtm_main_exch,
                    exchange_type="direct",
                    durable=True,
                    auto_delete=False,
                )
                channel.queue.declare(
                    queue=self.jtm_task_request_q,
                    durable=True,
                    exclusive=False,
                    auto_delete=True,
                )
                channel.queue.bind(
                    exchange=self.jgi_jtm_main_exch,
                    queue=self.jtm_task_request_q,
                    routing_key=self.jtm_task_request_q,
                )
                channel.basic.consume(
                    self.process_jtm_command, self.jtm_task_request_q, no_ack=False
                )
                channel.basic.qos(prefetch_count=1)
                channel.start_consuming()
                if not channel.consumer_tags:
                    channel.close()
            except amqpstorm.AMQPError as why:
                logger.exception(why)
                self.create_connection()
            except KeyboardInterrupt:
                self.connection.close()
                break
            except Exception as e:
                logger.exception(e)
                if self.connection:
                    self.connection.close()
                kill_child_proc(self.ppid)
                raise

    def process_jtm_command(self, message):
        """
        Callback for handling task request from JtmInterface

        :param message:
        :return:
        """
        msg_unzipped = json.loads(zloads(message.body))
        logger.info("New task: {}".format(msg_unzipped))
        task_type = msg_unzipped["task_type"]

        if task_type == "task":
            self.send_reply(message, "task", process_task_request(msg_unzipped))
        else:
            logger.critical("Task type not found: {}".format(task_type))
            self.send_reply(message, "Undefined task type", -1)


class JtmNonTaskCommandRunner(JtmAmqpstormBase):
    """
    Process non-task command request from JtmInterface <- JTM Client

    """

    def start(self):
        """Start the JtmNonTaskCommandRunner.
        :return:
        """
        if not self.connection:
            self.create_connection()
        while True:
            try:
                channel = self.connection.channel()
                channel.exchange.declare(
                    exchange=self.jgi_jtm_main_exch,
                    exchange_type="direct",
                    durable=True,
                    auto_delete=False,
                )
                channel.queue.declare(
                    queue=self.jtm_status_request_q,
                    durable=True,
                    exclusive=False,
                    auto_delete=True,
                )
                channel.queue.bind(
                    exchange=self.jgi_jtm_main_exch,
                    queue=self.jtm_status_request_q,
                    routing_key=self.jtm_status_request_q,
                )
                channel.basic.consume(
                    self.process_jtm_command, self.jtm_status_request_q, no_ack=False
                )
                channel.basic.qos(prefetch_count=1)
                channel.start_consuming()
                if not channel.consumer_tags:
                    channel.close()
            except amqpstorm.AMQPError as why:
                logger.exception(why)
                self.create_connection()
            except KeyboardInterrupt:
                self.connection.close()
                break
            except Exception as e:
                logger.exception(e)
                if self.connection:
                    self.connection.close()
                kill_child_proc(self.ppid)
                raise

    def process_jtm_command(self, message):
        """
        Callback for handling request from JtmInterface

        :param message:
        :return:
        """
        msg_unzipped = json.loads(zloads(message.body))
        logger.info("New task: {}".format(msg_unzipped))
        task_type = msg_unzipped["task_type"]

        # if task_type == "task":
        #     self.send_reply(message, "task", process_task_request(msg_unzipped))
        if task_type == "status":
            self.send_reply(
                message, "status", process_task_status(int(msg_unzipped["task_id"]))
            )
        elif task_type == "resource":
            self.send_reply(
                message, "resource", process_resource_log(int(msg_unzipped["task_id"]))
            )
        elif task_type == "kill":
            self.send_reply(
                message, "kill", process_task_kill(int(msg_unzipped["task_id"]))
            )
        elif task_type == "check_manager":
            self.send_reply(message, "check_manager", 88)
        elif task_type == "check_worker":
            self.send_reply(message, "check_worker", process_check_worker(msg_unzipped))
        elif task_type == "remove_pool":
            self.send_reply(
                message, "remove_pool", process_remove_pool(msg_unzipped["task_pool"])
            )
        else:
            logger.critical("Task type not found: {}".format(task_type))
            self.send_reply(message, "Undefined task type", -1)


# --------------------------------------------------------------------------------------------------
def extract_cromwell_run_id(task_id: int) -> str:
    """
    Extract Cromwell run id from "script" of a task

    :param task_id:
    :return: Cromwell run id in UUID format
    """
    db = DbSqlMysql(config=CONFIG)
    task_cmd_str = db.selectScalar(
        JTM_SQL["select_usercmd_tasks_by_taskid"] % dict(task_id=task_id), debug=False
    )
    db.close()
    # NOTE: here it is assumed that the script file full path exists at the
    # end of task_cmd_str
    try:
        task_script_file = task_cmd_str.split()[-1]
    except Exception as e:
        logger.exception("Invalid user cmd: %s" % e)
        task_script_file = None

    run_id = None
    # NOTE: here UUID format spec is assumed to comply with "8-4-4-4-12" format
    if os.path.isfile(task_script_file) and os.access(task_script_file, os.R_OK):
        with open(task_script_file, "r") as script_file:
            for line in script_file:
                run_id = extract_cromwell_id(line)
                if run_id:
                    break

    return run_id


# --------------------------------------------------------------------------------------------------
def send_update_task_status_msg(
    task_id: int, status_from, status_to: int, fail_code=None, reason=None
):
    """
    Publish a message for pushing task status change to JAWS Site

    :param task_id:
    :param status_from: None or status code 0 ~ 4 or -2 if failed
    :param status_to: status code 0 ~ 4 or -2 if failed
    :param fail_code: fail code if failed -1 ~ -7
    :param reason: additional info string like failure reason
    :return: None
    """
    now = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    run_id = extract_cromwell_run_id(task_id)
    reversed_task_status = dict(map(reversed, CONFIG.constants.TASK_STATUS.items()))
    reversed_done_flags = dict(map(reversed, CONFIG.constants.DONE_FLAGS.items()))
    reason_str = ""
    if fail_code:
        reason_str += "%s, " % reversed_done_flags[fail_code]
    if reason:
        reason_str += "%s" % reason

    data = {
        "cromwell_run_id": run_id,  # this is not the JAWS run_id
        "cromwell_job_id": task_id,
        "status_from": reversed_task_status[status_from]
        if status_from is not None
        else "",
        "status_to": reversed_task_status[status_to] if status_to is not None else "",
        "timestamp": now,
        "reason": reason_str,
    }

    # send message to Site
    try:
        with rpc_client.RpcClient(
            {
                "host": CONFIG.configparser.get("SITE_RPC_CLIENT", "host"),
                "vhost": CONFIG.configparser.get("SITE_RPC_CLIENT", "vhost"),
                "port": CONFIG.configparser.get("SITE_RPC_CLIENT", "port"),
                "user": CONFIG.configparser.get("SITE_RPC_CLIENT", "user"),
                "queue": CONFIG.configparser.get("SITE_RPC_CLIENT", "queue"),
                "password": CONFIG.configparser.get("SITE_RPC_CLIENT", "password"),
            }
        ) as rpc_cl:
            wait_count = 0
            response = rpc_cl.request("update_job_status", data)
            logger.debug(f"Return msg from JAWS Site: {response}")
            while (
                "error" in response and response["error"]["message"] == "Server timeout"
            ):
                wait_count += 1
                if wait_count == 60:  # try for 1min
                    logger.error("RPC reply timeout!")
                    break
                logger.debug(
                    f"RPC reply delay. Wait for a result from JAWS Site RPC server: {response}"
                )
                time.sleep(1.0)
                response = rpc_cl.request("update_job_status", data)
    except Exception as error:
        logger.error(f"RPC call failed: {error}")
        raise

    if "result" in response:
        logger.debug(f"Status change message sent successfully: {data}")
        pass
    else:
        logger.error(f"Status update failed: {response['error']['message']}")


# --------------------------------------------------------------------------------------------------
def recv_hb_from_worker_proc(hb_queue_name, log_dest_dir, b_resource_log):
    """

    :param hb_queue_name:
    :param log_dest_dir:
    :param b_resource_log:
    :return:
    """
    hb_msg = CONFIG.constants.HB_MSG
    jtm_worker_hb_exch = CONFIG.configparser.get("JTM", "jtm_worker_hb_exch")
    interval = CONFIG.configparser.getfloat("JTM", "client_hb_recv_interval")
    hb_mgs_to_recv = CONFIG.configparser.getint("JTM", "heartbeat_message_count")

    with RmqConnectionAmqpstorm(config=CONFIG).open() as conn:
        with conn.channel() as ch:
            ch.queue.declare(
                queue=hb_queue_name, durable=False, exclusive=False, auto_delete=True
            )
            ch.exchange.declare(
                exchange=jtm_worker_hb_exch,
                exchange_type="direct",
                durable=False,
                auto_delete=False,
            )
            ch.queue.bind(
                exchange=jtm_worker_hb_exch,
                queue=hb_queue_name,
                routing_key=hb_queue_name,
            )

            b_is_msg_cleared = False  # all stacked messages are processed or not
            b_is_worker_found = False  # is any alive worker
            max_worker_check_count = 0  # max number of checking workers

            while True:
                worker_ids_dict = {}
                # this qos is to set how many messages to take out from RMQ queue
                ch.basic.qos(hb_mgs_to_recv)
                message = ch.basic.get(queue=hb_queue_name, no_ack=True)
                if message and not b_is_msg_cleared:
                    cnt = message.method["message_count"]
                    logger.debug(f"hb msg count = {cnt}")
                    for i in range(cnt):
                        _ = ch.basic.get(queue=hb_queue_name, no_ack=True)
                    b_is_msg_cleared = True

                elif message and b_is_msg_cleared:
                    cnt = message.method["message_count"]
                    logger.debug(f"hb msg count = {cnt}")
                    msg_unzipped = json.loads(zloads(message.body))
                    msg_unzipped = {int(k): v for k, v in msg_unzipped.items()}
                    a_worker_id = msg_unzipped[hb_msg["worker_id"]]
                    worker_ids_dict[a_worker_id] = msg_unzipped

                    for i in range(cnt):
                        message = ch.basic.get(queue=hb_queue_name, no_ack=True)
                        if message:
                            msg_unzipped = json.loads(zloads(message.body))
                            msg_unzipped = {int(k): v for k, v in msg_unzipped.items()}
                            a_worker_id = msg_unzipped[hb_msg["worker_id"]]
                            worker_ids_dict[a_worker_id] = msg_unzipped

                    for k, v in worker_ids_dict.items():
                        task_id = v[hb_msg["task_id"]]
                        root_proc_id = v[hb_msg["root_pid"]]
                        child_proc_id = v[hb_msg["child_pid"]]
                        a_worker_id = v[hb_msg["worker_id"]]
                        slurm_job_id = v[hb_msg["slurm_jobid"]]
                        worker_type = v[hb_msg["worker_type"]]
                        end_datetime = v[hb_msg["end_date"]]
                        life_left = v[hb_msg["life_left"]]
                        mem_per_node = v[hb_msg["mem_per_node"]]
                        mem_per_core = v[hb_msg["mem_per_core"]]
                        num_cores = v[hb_msg["num_cores"]]
                        job_time = v[hb_msg["job_time"]]
                        clone_time = v[hb_msg["clone_time_rate"]]
                        host_name = v[hb_msg["host_name"]]
                        jtm_host_name = v[hb_msg["jtm_host_name"]]
                        ip_addr = v[hb_msg["ip_address"]]
                        pool_name = v[hb_msg["pool_name"]]

                        if b_resource_log:
                            logger.resource(v)
                        if pool_name:
                            queue_name = (
                                CONFIG.configparser.get("JTM", "jtm_inner_request_q")
                                + "."
                                + pool_name
                            )
                        else:
                            queue_name = ""

                        success = False

                        # Todo: still need to do "while" for checking db connection?
                        while success is not True:
                            success = True
                            try:
                                db = DbSqlMysql(config=CONFIG)
                                bExists = db.selectScalar(
                                    JTM_SQL["select_exists_workers_by_workerid"]
                                    % dict(worker_id=a_worker_id)
                                )

                                if not bExists:
                                    # Todo: remove num_worker_on_the_node insertion
                                    db.execute(
                                        JTM_SQL["insert_workers_workerid_slurmjobid"]
                                        % dict(
                                            worker_id=a_worker_id,
                                            slurm_jid=slurm_job_id,
                                            worker_type=worker_type,
                                            nwpn=1,
                                        )
                                    )
                                else:
                                    # Todo: stress test needed! Test if nworkers > 1000
                                    # Todo: can only update end_datetime and life_left after 1st insert
                                    db.execute(
                                        JTM_SQL[
                                            "update_workers_enddate_lifeleft_by_workerid"
                                        ]
                                        % dict(
                                            worker_id=a_worker_id,
                                            now=end_datetime,
                                            life_left=life_left,
                                            mpn=mem_per_node if worker_type != 0 else 0,
                                            mpc=mem_per_core if worker_type != 0 else 0,
                                            num_cores=num_cores
                                            if worker_type != 0
                                            else 0,
                                            job_time=job_time,
                                            clone_rate=clone_time
                                            if worker_type != 0
                                            else 0,
                                            host_name=host_name,
                                            jtm_host_name=jtm_host_name,
                                            ipaddr=ip_addr,
                                            pool_name=queue_name,
                                            slurm_jid=slurm_job_id,
                                        )
                                    )

                                db.commit()
                                db.close()
                            except Exception as e:
                                logger.critical(e)
                                logger.critical(
                                    "Failed to update workers table for enddate and lifeleft."
                                )
                                logger.debug(
                                    "Retry to update workers table for enddate and lifeleft."
                                )
                                success = False
                                conn.sleep(1)

                        # handle "pending"
                        if (
                            task_id > 0
                            and root_proc_id == child_proc_id
                            and slurm_job_id > 0
                        ):
                            pass
                        elif (
                            task_id > 0 and root_proc_id != child_proc_id
                        ):  # if user process processing started
                            logger.debug("Task %d status ==> running" % task_id)
                            datetime_str = datetime.datetime.now().strftime("%Y-%m-%d")
                            if log_dest_dir:
                                log_dir_name = os.path.join(log_dest_dir, "resource")
                            else:
                                log_dir_name = "%s/resource" % (os.getcwd())

                            padded_dir_str = pad_string_path(
                                task_id, depth=3
                            )  # task id 226 --> 00/00/02
                            make_dir(os.path.join(log_dir_name, padded_dir_str))
                            resource_log_fname = "%s/%s/jtm_resource_%d_%s.log" % (
                                log_dir_name,
                                padded_dir_str,
                                task_id,
                                datetime_str,
                            )

                            with open(resource_log_fname, "a") as rf:
                                rf.write(",".join([str(i) for i in v.values()]))
                                rf.write("\n")

                            os.chmod(resource_log_fname, 0o777)

                            try:
                                # Update tasks table with "running" status == 2 if status is still 0 or 1
                                db = DbSqlMysql(config=CONFIG)
                                status_now = db.selectScalar(
                                    JTM_SQL["select_status_runs_by_taskid"]
                                    % dict(task_id=task_id),
                                    debug=False,
                                )
                                logger.debug(f"status now {task_id} = {status_now}")

                                if not STANDALONE:
                                    logger.debug(
                                        f"status change msg {task_id}: {status_now} => running"
                                    )
                                    if status_now == TASK_STATUS["pending"]:
                                        send_update_task_status_msg(
                                            task_id,
                                            status_now,
                                            TASK_STATUS["running"],
                                            reason=slurm_job_id,
                                        )

                                db.execute(
                                    JTM_SQL["update_runs_status_to_running_by_taskid"]
                                    % dict(
                                        status_id=TASK_STATUS["running"],
                                        task_id=task_id,
                                        worker_id=a_worker_id,
                                        child_pid=child_proc_id,
                                        resource_log=resource_log_fname,
                                    ),
                                    debug=False,
                                )
                                db.commit()
                                db.close()
                            except Exception as e:
                                logger.critical(e)
                                logger.critical(
                                    "Failed to update runs table for status."
                                )
                                ch.close()
                                conn.close()
                                raise

                    # Collect worker_ids from hb
                    alive_worker_id_list = []
                    for k, v in worker_ids_dict.items():
                        alive_worker_id_list.append(v[hb_msg["worker_id"]])

                    # Collect worker_id which are still set as alive from workers table
                    success = False
                    while success is not True:
                        success = True
                        try:
                            db = DbSqlMysql(config=CONFIG)
                            live_worker_id_list = db.selectAll(
                                JTM_SQL[
                                    "select_workerid_workers_by_lifeleft_jtmhostname"
                                ]
                                % dict(jtm_host_name=jtm_host_name),
                                debug=False,
                            )
                            live_worker_id_list = [i[0] for i in live_worker_id_list]

                            # Check & update workers table
                            # Set life_left as -1 for dead workers
                            for w in live_worker_id_list:
                                if w not in alive_worker_id_list:
                                    # Update lifeleft --> -1
                                    db.execute(
                                        JTM_SQL["update_workers_lifeleft_by_workerid"]
                                        % dict(worker_id=w, jtm_host_name=jtm_host_name)
                                    )
                            db.commit()
                            db.close()
                        except Exception as e:
                            logger.critical(e)
                            logger.critical(
                                "Failed to update workers table for lifeleft."
                            )
                            logger.debug("Retry to update workers table for lifeleft.")
                            success = False
                            conn.sleep(1)

                    db = DbSqlMysql(config=CONFIG)
                    alive_total_num_workers = db.selectScalar(
                        JTM_SQL["select_sum_nwpn_workers_by_lifeleftt_jtmhostname"]
                        % dict(jtm_host_name=jtm_host_name),
                        debug=False,
                    )
                    db.close()
                    NUM_TOTAL_WORKERS.value = (
                        int(alive_total_num_workers) if alive_total_num_workers else 0
                    )
                    logger.debug(
                        "# workers: in hb=%d, in table=%d, alive+requested=NUM_TOTAL_WORKERS=%d"
                        % (
                            len(alive_worker_id_list),
                            len(live_worker_id_list),
                            NUM_TOTAL_WORKERS.value,
                        )
                    )

                    if NUM_TOTAL_WORKERS.value > 0:
                        b_is_worker_found = True
                        max_worker_check_count = 0  # reinitialize

                elif not message and b_is_msg_cleared and b_is_worker_found:
                    # If we are here, unfortunately, we lost all the workers that we've been using.
                    max_worker_check_count += 1
                    NUM_TOTAL_WORKERS.value = 0
                    logger.info("Waiting for worker(s)...")
                    hb_check = CONFIG.configparser.getint(
                        "JTM", "worker_hb_check_max_count"
                    )
                    if (
                        hb_check != 0 and max_worker_check_count > hb_check
                    ):  # hit the max checking limit
                        # Close connection and kill parent and itself
                        ch.close()
                        conn.close()
                        raise OSError(2, "Worker not found")

                    # If there no workers alive after 5 checks, set life_left to -1 for all
                    if max_worker_check_count == 3:
                        try:
                            db = DbSqlMysql(config=CONFIG)
                            db.execute(JTM_SQL["update_workers_lifeleft_for_last"])
                            db.commit()
                            db.close()
                        except Exception as e:
                            logger.critical(e)
                            logger.critical(
                                "Failed to update workers table for lifeleft."
                            )
                            ch.close()
                            conn.close()
                            raise

                time.sleep(interval)


# -------------------------------------------------------------------------------
def process_task_request(msg):
    """
    Get task request from jtm-submit and send it to a worker

    :param msg: unzipped dict
    :param inner_task_request_queue: inner task queue name (jgi-task-manager --> worker)
    :return:
    """
    # Parse the msg from jtm-submit
    if "command" not in msg:
        logger.critical("Critical: cannot find user command.")
        return -5

    if "task_type" not in msg:
        logger.critical("Critical: cannot find task type.")
        return -5

    user_task_cmd = msg["command"]
    task_type = msg["task_type"]
    output_file = (
        msg["output_files"] if "output_files" in msg else ""
    )  # comma separated list ex) "a.out,b.out,c.out"
    output_dir = msg.get("output_dir", "")
    stdout_file = msg.get("stdout", "")
    stderr_file = msg.get("stderr", "")

    # If pool is set, the tasks which use the queue name (=pool name)
    # will only be sent to the pool of workers
    jtm_inner_request_q = CONFIG.configparser.get("JTM", "jtm_inner_request_q")
    inner_task_request_queue = jtm_inner_request_q

    if "pool" in msg:
        pool_spec_json_str = json.loads(json.dumps(msg["pool"]))
        pool_name = None
        if "name" in pool_spec_json_str and pool_spec_json_str["name"]:
            pool_name = pool_spec_json_str["name"]
        if pool_name:
            # _jtm_inner_request_queue.<cluster_name>.jtm.<pool_name>
            inner_task_request_queue = jtm_inner_request_q + "." + pool_name
        else:
            inner_task_request_queue = jtm_inner_request_q + ".small"

    # Deal with custom pool of workers
    last_task_id = -1
    b_failed_to_request_worker = False
    w_int = CONFIG.configparser.getfloat("JTM", "worker_hb_recv_interval")
    slurm_job_id = 0

    if "pool" in msg and "name" in msg["pool"] and "time" in msg["pool"]:
        # {u'resource': u'cori', u'name': u'test', u'size': 1}
        # ==> {"resource": "cori", "name": "test", "size": 1}
        pool_spec_json_str = json.loads(json.dumps(msg["pool"]))

        # Set default values defined in Config.py
        assert "name" in pool_spec_json_str
        pool_name = pool_spec_json_str["name"]
        assert pool_name is not None
        pool_cluster = CONFIG.configparser.get("JTM", "cluster")
        pool_ncpus = CONFIG.configparser.getint("SLURM", "ncpus")
        pool_mem = CONFIG.configparser.get("SLURM", "mempernode")

        # if not defined in the configuration file, the below will have empty string value
        pool_constraint = CONFIG.configparser.get("SLURM", "constraint")
        pool_charge_account = CONFIG.configparser.get("SLURM", "charge_accnt")
        pool_qos = CONFIG.configparser.get("SLURM", "qos")
        pool_partition = CONFIG.configparser.get("SLURM", "partition")

        # Note: pool size = num_nodes_to_request * num_workers_per_node
        num_workers_per_node = CONFIG.configparser.getint("JTM", "num_workers_per_node")
        num_nodes_to_request = CONFIG.configparser.getint("SLURM", "nnodes")

        # Worker type is restricted to "dynamic" for now.
        # Todo: add the feature to remove the custom pool by user task to support "static"
        #   workers in custom pool
        # NOTE: if pool_size is set, which means user sets the number of workers needed,
        # n number of jtm dynamic workers will be created, and the pool of n dynamic
        # workers won't be increased over n. The workers will be terminated if there is
        # no tasks requested
        # if not,
        # only one dynamic worker will be created. The worker will be terminated if there is
        # no tasks in the queue for a specified time duration
        pool_cluster = pool_spec_json_str.get("cluster")
        pool_time = pool_spec_json_str.get("time")
        pool_ncpus = pool_spec_json_str.get("cpu")
        pool_mem = pool_spec_json_str.get("mem")
        pool_constraint = pool_spec_json_str.get("constraint", "")
        pool_qos = pool_spec_json_str.get("qos", "")
        pool_partition = pool_spec_json_str.get("partition", "")
        pool_charge_account = pool_spec_json_str.get("account", "")
        num_workers_per_node = pool_spec_json_str.get("nwpn", 1)
        num_nodes_to_request = pool_spec_json_str.get("node", 1)

        assert len(user_task_cmd) <= 1024
        assert len(output_file) <= 1024

        ################
        # Create a pool
        ################
        # Steps
        # 1. check the number if workers in the custom pool
        # 2. if 0, create a pool
        #    else send tasks to the pool
        db = DbSqlMysql(config=CONFIG)
        num_slurm_jid_in_pool = db.selectScalar(
            JTM_SQL["select_count_distinct_jid_workers_by_poolname"]
            % dict(pool_name=inner_task_request_queue, hbinterval=w_int * 3),
            debug=False,
        )

        # Note: slurm jid --> sacct --> doublecheck the status of node allocation and worker
        #  if no real alive workers --> request new node
        #  if num_worker_to_add == 0, squeue -j slurm jid
        ####################################################################################
        slurm_jid_list = db.selectAll(
            JTM_SQL["select_distinct_jid_workers_by_poolname"]
            % dict(pool_name=inner_task_request_queue, hbinterval=w_int * 3)
        )
        db.close()

        slurm_jid_list = [i[0] for i in slurm_jid_list]
        logger.debug(
            "slurm job id for {}: {}".format(inner_task_request_queue, slurm_jid_list)
        )

        for jid in slurm_jid_list:
            cmd = "squeue -j %d" % jid
            so, _, ec = run_sh_command(cmd, log=logger, show_stdout=False)
            if not re.search(r"\bPD\b", so) and not re.search(r"\bR\b", so):
                num_slurm_jid_in_pool -= 1

                # Delete worker info
                logger.info("Found dead worker(s) from %s, %s" % (str(jid), so))
                db = DbSqlMysql(config=CONFIG)
                db.execute(
                    JTM_SQL["delete_from_workers_by_slurmjid"]
                    % dict(
                        slurm_jid=jid,
                    ),
                    debug=DEBUG,
                )
                db.commit()
                db.close()
        ####################################################################################

        pool_size = num_nodes_to_request * num_workers_per_node

        # This is the actual number of workers = jid * nwpn
        num_live_worker_in_pool = num_slurm_jid_in_pool * num_workers_per_node
        num_worker_to_add = int(
            ceil(float(pool_size - num_live_worker_in_pool) / num_workers_per_node)
        )

        logger.debug(
            "num_worker_to_add={} pool_size={} "
            "num_live_worker_in_pool={} num_slurm_jid_in_pool={} "
            "NUM_TOTAL_WORKERS={} num_workers_per_node={}"
            "".format(
                num_worker_to_add,
                pool_size,
                num_live_worker_in_pool,
                num_slurm_jid_in_pool,
                NUM_TOTAL_WORKERS.value,
                num_workers_per_node,
            )
        )

        uniq_worker_id = None
        env_act = CONFIG.configparser.get("JTM", "env_activation")
        for i in range(0, num_worker_to_add):
            b_failed_to_request_worker = False
            uniq_worker_id = str(shortuuid.uuid())
            sbatch_cmd_str = """{}jtm {} worker \
                -wt dynamic \
                -p {} \
                -cl {} \
                -c {} \
                -t {} \
                -m {} \
                -wi {} \
                -nwpn {} \
                {} {} {} {}""".format(
                "%s && " % env_act if env_act else "",
                "--config=%s" % CONFIG.config_file if CONFIG else "",
                pool_name,
                pool_cluster,
                pool_ncpus,
                pool_time,
                pool_mem,
                uniq_worker_id,
                num_workers_per_node,
                "-C %s" % pool_constraint if pool_constraint else "",
                "--qos %s" % pool_qos if pool_qos else "",
                "-A %s" % pool_charge_account if pool_charge_account else "",
                "-P %s" % pool_partition if pool_partition else "",
            )

            logger.info("Executing sbatch {}".format(sbatch_cmd_str))

            # Run sbatch from jtm worker command
            so, _, ec = run_sh_command(sbatch_cmd_str, log=logger)

            # Print batch job script for logging
            run_sh_command(sbatch_cmd_str + " --dry_run", log=logger)

            # Get the slurm job id returned from jtm-worker
            try:
                regex = re.compile(r"(\d+)", re.I)
                match = regex.search(so)
                slurm_job_id = int(match.group(1))
            except Exception:
                logger.critical(
                    "Failed to get a valid job ID back from requesting a dynamic worker"
                )
                ec = 1  # make it fail

            if ec == 0:  # if sbatch by jtm-worker done successfully
                for serial_worker_num in range(num_workers_per_node):
                    # Insert into workers for new worker id
                    try:
                        db = DbSqlMysql(config=CONFIG)
                        # If a dynamic worker is requested successfully,
                        # insert the info into workers table
                        # logger.info("Try to update workers table")check_worker

                        # Fixme: After sbatch, another sbatch with the same pool name will be executed
                        #  again.
                        # Solution: Set the lifeleft=-2 and update select_count_workers_by_poolname sql
                        # statement to check only lifeleft!=-1 so that num_live_worker_in_pool can include
                        # already sbatched workers for the pool.
                        #
                        db.execute(
                            JTM_SQL["insert_workers_workerid_workertype_poolname"]
                            % dict(
                                worker_id=uniq_worker_id + str(serial_worker_num + 1),
                                worker_type=CONFIG.constants.WORKER_TYPE["dynamic"],
                                pool_name=inner_task_request_queue,
                                jtm_host_name=pool_cluster,
                                lifeleft=-2,
                                slurm_jid=slurm_job_id,
                                nwpn=1,
                            )
                        )
                        db.commit()
                        db.close()
                    except Exception as e:
                        logger.critical(
                            "Failed to insert workers table for workerid and workertype."
                        )
                        logger.critical(e)
                        logger.debug(
                            "Retry to insert workers table for workerid and workertype."
                        )
                        raise
            else:
                logger.critical("Failed to execute the command, %s" % (sbatch_cmd_str))
                b_failed_to_request_worker = True
                last_task_id = TASK_STATUS["invalidtask"]
                break

    if not b_failed_to_request_worker:
        last_task_id = -1

        # Fixme: mysql connection.py IndexError: bytearray index out of range
        # note: seems like mysql pool connection error
        success = False
        while success is not True:
            success = True
            try:
                db = DbSqlMysql(config=CONFIG)
                # table fields: userCmd, outFiles, doneFlag, retryCnt, task_type
                db.execute(
                    JTM_SQL["insert_tasks_usercmd_outfiles"]
                    % (
                        user_task_cmd,
                        output_file,
                        "0",
                        0,
                        CONFIG.constants.TASK_TYPE[task_type],
                    )
                )
                db.commit()
                last_task_id = db.selectScalar(JTM_SQL["select_last_insert_id"])
                db.close()
                logger.debug("last_task_id = %d" % last_task_id)
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed to insert tasks table for new user task.")
                logger.debug("Retry to insert tasks table for a user task.")
                success = False
                time.sleep(1.0)

        success = False
        while success is not True:
            success = True
            try:
                logger.debug("Task %d status ==> ready" % last_task_id)
                db = DbSqlMysql(config=CONFIG)
                db.execute(
                    JTM_SQL["insert_runs_tid_sid"]
                    % dict(task_id=last_task_id, status_id=TASK_STATUS["ready"])
                )
                db.commit()
                db.close()
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed to insert runs table for a new run.")
                logger.debug("Retry to insert runs table for a new run.")
                success = False
                time.sleep(1.0)

        # if it successfully updates runs table and gets a task id
        if last_task_id != -1:
            db = DbSqlMysql(config=CONFIG)
            if not STANDALONE:
                logger.debug("process task request -->")
                logger.debug(f"status change msg {last_task_id}: created => ready")
                send_update_task_status_msg(
                    last_task_id, TASK_STATUS["created"], TASK_STATUS["ready"]
                )

            task_status_int = int(
                db.selectScalar(
                    JTM_SQL["select_status_runs_by_taskid"] % dict(task_id=last_task_id)
                )
            )
            b_already_canceled = int(
                db.selectScalar(
                    JTM_SQL["select_cancelled_runs_by_taskid"]
                    % dict(task_id=last_task_id)
                )
            )
            db.close()
            logger.debug(f"current task status = {task_status_int}")

            # Note: this is just in case
            #  It is based on the assumption that a task can be cancelled between ready -> queued
            #  status change.
            status_now = None
            if b_already_canceled == 1 or task_status_int == TASK_STATUS["terminated"]:
                if task_status_int != TASK_STATUS["terminated"]:
                    db = DbSqlMysql(config=CONFIG)
                    db.execute(
                        JTM_SQL["update_tasks_doneflag_by_taskid"]
                        % dict(
                            task_id=last_task_id, done_flag=TASK_STATUS["terminated"]
                        )
                    )
                    db.commit()
                    db.close()

                    if not STANDALONE:
                        logger.debug("process task request -->")
                        logger.debug(
                            f"status change msg {last_task_id}: {task_status_int} => terminated"
                        )
                        send_update_task_status_msg(
                            last_task_id,
                            task_status_int,
                            TASK_STATUS["failed"],
                            fail_code=TASK_STATUS["terminated"],
                        )
            else:
                # Prepare msg to jtm-worker
                msg_to_send_dict = {}
                msg_to_send_dict["task_id"] = last_task_id
                msg_to_send_dict["user_cmd"] = user_task_cmd
                msg_to_send_dict["output_files"] = output_file
                msg_to_send_dict["done_flag"] = 0
                msg_to_send_dict["task_type"] = CONFIG.constants.TASK_TYPE[task_type]
                msg_to_send_dict["output_dir"] = output_dir
                msg_to_send_dict["stdout"] = stdout_file
                msg_to_send_dict["stderr"] = stderr_file
                logger.info(
                    "Total number of workers (alive + requested): %d",
                    NUM_TOTAL_WORKERS.value,
                )

                # Create and send request message to workers
                msg_zipped = zdumps(json.dumps(msg_to_send_dict))
                corr_id = str(uuid.uuid4())
                exch_name = CONFIG.configparser.get("JTM", "jtm_inner_main_exch")
                reply_q = CONFIG.configparser.get("JTM", "jtm_inner_result_q")
                logger.info("Send a task to {}".format(inner_task_request_queue))

                # Just a little gap b/w ready and queued
                # time.sleep(CONFIG.configparser.getfloat("JTM", "task_stat_update_interval"))

                # This TASK_STATUS["queued"] means it's queued to jtm task queue
                # TASK_STATUS["pending"] means node requested to slurm
                # TASK_STATUS["running"] means task processing started
                try:
                    logger.debug("Task %d status ==> queued" % last_task_id)
                    db = DbSqlMysql(config=CONFIG)
                    status_now = db.selectScalar(
                        JTM_SQL["select_status_runs_by_taskid"]
                        % dict(task_id=last_task_id),
                        debug=False,
                    )
                    logger.debug(f"status now {last_task_id} = {status_now}")

                    if not STANDALONE:
                        logger.debug("process task request -->")
                        logger.debug(
                            f"status change msg {last_task_id}: {status_now} => queued"
                        )
                        if (
                            status_now in (TASK_STATUS["ready"], TASK_STATUS["success"])
                            or status_now < 0
                        ):
                            send_update_task_status_msg(
                                last_task_id, status_now, TASK_STATUS["queued"]
                            )

                    # Note: only when the previous status == ready
                    if status_now == TASK_STATUS["ready"]:
                        db.execute(
                            JTM_SQL["update_runs_status_by_taskid"]
                            % dict(
                                status_id=TASK_STATUS["queued"],
                                now=time.strftime("%Y-%m-%d %H:%M:%S"),
                                task_id=last_task_id,
                            ),
                            debug=False,
                        )
                        db.commit()

                    status_now = db.selectScalar(
                        JTM_SQL["select_status_runs_by_taskid"]
                        % dict(task_id=last_task_id),
                        debug=False,
                    )
                    logger.debug(f"status now {last_task_id} = {status_now}")
                    db.close()
                except Exception as e:
                    logger.critical(e)
                    logger.critical(
                        "Failed to update runs table for status and startdate."
                    )
                    # Todo: properly update runs for the failure
                    # Todo: set the task status --> failed
                    raise

                # now waiting for slurm allocation
                db = DbSqlMysql(config=CONFIG)
                status_now = db.selectScalar(
                    JTM_SQL["select_status_runs_by_taskid"]
                    % dict(task_id=last_task_id),
                    debug=False,
                )
                logger.debug(f"status now {last_task_id} = {status_now}")

                if not STANDALONE:
                    logger.debug("process task request -->")
                    logger.debug(
                        f"status change msg {last_task_id}: {status_now} => pending"
                    )
                    if (
                        status_now in (TASK_STATUS["queued"], TASK_STATUS["success"])
                        or status_now < 0
                    ):
                        send_update_task_status_msg(
                            last_task_id,
                            status_now,
                            TASK_STATUS["pending"],
                            reason=slurm_job_id,
                        )

                # Note: only when the previous status == queued
                if status_now == TASK_STATUS["queued"]:
                    db.execute(
                        JTM_SQL["update_runs_status_by_taskid"]
                        % dict(
                            status_id=TASK_STATUS["pending"],
                            now=time.strftime("%Y-%m-%d %H:%M:%S"),
                            task_id=last_task_id,
                        ),
                        debug=False,
                    )
                    db.commit()

                status_now = db.selectScalar(
                    JTM_SQL["select_status_runs_by_taskid"]
                    % dict(task_id=last_task_id),
                    debug=False,
                )
                logger.debug(f"status now {last_task_id} = {status_now}")
                db.close()

                try:
                    with RmqConnectionAmqpstorm(config=CONFIG).open() as conn:
                        with conn.channel() as ch:
                            ch.exchange.declare(
                                exchange=exch_name,
                                exchange_type="direct",
                                durable=True,
                                auto_delete=False,
                            )
                            ch.queue.declare(
                                queue=inner_task_request_queue,
                                durable=True,
                                exclusive=False,
                                auto_delete=True,
                            )
                            ch.queue.bind(
                                exchange=exch_name,
                                queue=inner_task_request_queue,
                                routing_key=inner_task_request_queue,
                            )
                            properties = {
                                "correlation_id": corr_id,
                                "reply_to": reply_q,
                            }
                            response = amqpstorm.Message.create(
                                ch, msg_zipped, properties
                            )
                            response.publish(inner_task_request_queue)
                except Exception as detail:
                    logger.exception(
                        "Exception: Failed to send a request to %s",
                        inner_task_request_queue,
                    )
                    logger.exception("Detail: %s", str(detail))
                    raise OSError(2, "Failed to send a request to a worker")

        return last_task_id


# -------------------------------------------------------------------------------
def process_task_status(task_id):
    """
    Process task status request

    :param task_id: for task ID
    :return:
    """
    task_status_int = 0
    db = DbSqlMysql(config=CONFIG)
    cur = db.execute(
        JTM_SQL["select_status_cancelled_runs_by_taskid"] % dict(task_id=task_id)
    )
    ret = cur.fetchone()
    db.close()

    if ret is not None:
        logger.debug(
            "Task status from runs table: status, cancelled = {}, "
            "Task id = {}"
            "".format(ret, task_id)
        )
        if ret[1] == 1:  # cancellation requested
            task_status_int = TASK_STATUS["terminated"]
        else:
            # 0 if ready
            # 1 if queued
            # 2 if pending
            # 3 if running
            task_status_int = ret[0]
    else:
        logger.debug("Task ID not found: {}".format(task_id))
        task_status_int = TASK_STATUS["invalidtask"]  # invalid task id

    return task_status_int


# -------------------------------------------------------------------------------
def process_resource_log(task_id):
    """
    With given task id, get the resource log file and send to jtm-resource-log

    :param task_id:
    :return:
    """
    resource_log_file = None
    db = DbSqlMysql(config=CONFIG)
    cur = db.execute(JTM_SQL["select_resource_runs_by_taskid"] % dict(task_id=task_id))
    ret = cur.fetchone()
    db.close()
    if ret:
        resource_log_file = ret[0]

    return resource_log_file


# -------------------------------------------------------------------------------
def send_task_kill_request(task_id, wid, cpid):
    """
    Send a task kill request
    1. if wid is ready, send wid with task_id
    2. else
           any worker is not ready, so just send task_id

    :param task_id: task id
    :param wid: worker id
    :param cpid: child process id
    :return:
    """
    exch_name = CONFIG.configparser.get("JTM", "jtm_task_kill_exch")
    queue_name = CONFIG.configparser.get("JTM", "jtm_task_kill_q")
    message = {"task_id": task_id, "worker_id": wid, "child_pid": cpid}
    msg_zipped = zdumps(json.dumps(message))

    try:
        with RmqConnectionAmqpstorm(config=CONFIG).open() as conn:
            with conn.channel() as ch:
                ch.exchange.declare(
                    exchange=exch_name,
                    exchange_type="fanout",
                    durable=True,
                    auto_delete=False,
                )
                ch.queue.declare(
                    queue=queue_name, durable=True, exclusive=False, auto_delete=True
                )
                ch.queue.bind(exchange=exch_name, queue=queue_name, routing_key=wid)
                # Don't need corr_if and reply_q b/c the consumer won't
                # send a reply back to this
                properties = {}
                response = amqpstorm.Message.create(ch, msg_zipped, properties)
                response.publish(queue_name)
    except Exception as detail:
        logger.exception("Exception: Failed to send a kill to %s", queue_name)
        logger.exception("Detail: %s", str(detail))
        raise OSError(2, "Failed to send a request to a worker")


# -------------------------------------------------------------------------------
def process_task_kill(task_id):
    """
    Process request to terminate a task

    Steps
    # Parse the msg from jtm-kill
    # Get workerid
    # Delete tasks row for task_id (runs table will be cascaded)
    # Send kill msg

    :param: task_id
    :return:
    """

    # task_id = int(msg["task_id"])
    return_val = None

    db = DbSqlMysql(config=CONFIG)
    task_status = db.selectScalar(
        JTM_SQL["select_status_runs_by_taskid"] % dict(task_id=task_id)
    )
    db.close()

    def update_runs_cancelled(task_id, status_now):
        db = DbSqlMysql(config=CONFIG)
        logger.debug(
            db.selectAll(
                "select taskId, status, cancelled from runs where taskId = %(task_id)s"
                % dict(
                    task_id=task_id,
                )
            )
        )

        db.execute(
            JTM_SQL["update_runs_cancelled_by_tid"]
            % dict(task_id=task_id, now=time.strftime("%Y-%m-%d %H:%M:%S")),
            debug=DEBUG,
        )
        db.commit()
        logger.debug(
            db.selectAll(
                "select taskId, status, cancelled from runs where taskId = %(task_id)s"
                % dict(
                    task_id=task_id,
                )
            )
        )
        logger.debug(
            "taskId, workerId2, childPid, cancelled = {}".format(
                db.selectAll(
                    "select taskId, workerId2, childPid, cancelled "
                    "from runs "
                    "where cancelled = 1 AND workerId2 != 0 AND childPid != 0"
                )
            )
        )
        db.close()
        # Note: why can't just set the "status" field to -4 = terminated here?
        # because worker id and pid is not known yet

        if not STANDALONE:
            logger.debug("task kill -->")
            logger.debug(f"status change msg {task_id}: {status_now} => terminated")
            send_update_task_status_msg(
                task_id,
                status_now,
                TASK_STATUS["failed"],
                fail_code=TASK_STATUS["terminated"],
            )

    if task_status:
        if task_status in (TASK_STATUS["ready"], TASK_STATUS["queued"]):
            logger.debug(
                "Task cancellation requested but the task, %d has already been queued."
                % (task_id)
            )
            logger.debug("The task will be terminated once it's started.")
            update_runs_cancelled(task_id, task_status)
            return_val = 0
        elif task_status == TASK_STATUS["pending"]:
            update_runs_cancelled(task_id, task_status)
            return_val = 0
        elif task_status == TASK_STATUS["running"]:
            logger.debug(
                "Task cancellation requested. The task, %s is being terminated."
                % (task_id)
            )
            update_runs_cancelled(task_id, task_status)

            # Get the wid and child pid
            db = DbSqlMysql(config=CONFIG)
            worker_id_list = db.selectAs1Col(
                JTM_SQL["select_workerid_workers_by_tid"] % dict(task_id=task_id),
                debug=False,
            )
            logger.debug("worker list = {}".format(worker_id_list))
            child_proc_id = db.selectScalar(
                JTM_SQL["select_chilpid_runs_by_tid"] % dict(task_id=task_id),
                debug=False,
            )
            logger.debug("child_proc_id = {}".format(child_proc_id))
            db.close()

            if child_proc_id > 0 and worker_id_list is not None:
                # Send task id and process id to worker id
                logger.info(
                    "Send task kill command: {} {} {}".format(
                        task_id, worker_id_list, child_proc_id
                    )
                )
                for wid in worker_id_list:
                    send_task_kill_request(task_id, wid, child_proc_id)

            return_val = 0
        elif task_status == TASK_STATUS["success"]:
            logger.debug(
                "Task cancellation requested but the task is in completed status."
            )
            return_val = 0
        elif task_status == TASK_STATUS["terminated"]:
            logger.debug(
                "Task cancellation requested but the task is already in terminated status."
            )
            return_val = 0
        elif task_status == TASK_STATUS["failed"]:
            logger.debug(
                "Task cancellation requested but the task is already in failed status."
            )
            return_val = 0
        elif task_status in (
            TASK_STATUS["outputerror"],
            TASK_STATUS["outputerror"],
            TASK_STATUS["outputerror"],
        ):
            logger.debug("Task cancellation request is ignored.")
            return_val = 0
        else:
            logger.debug("Failed to cancel a task. Unexpected condition.")
            return_val = -1
    else:  # task id not found
        return_val = -5

    return return_val


# -------------------------------------------------------------------------------
def process_check_worker(msg_unzipped):
    """

    :param msg_unzipped:
    :return:
    """
    # If custom_pool name is set, get the number of workers in the pool
    # else the total number of live workers will be sent
    w_int = CONFIG.configparser.getfloat("JTM", "worker_hb_recv_interval")
    jtm_inner_request_q = CONFIG.configparser.get("JTM", "jtm_inner_request_q")
    ret_val = 0

    if (
        "task_pool" in msg_unzipped
        and msg_unzipped["task_pool"]
        and "jtm_host_name" in msg_unzipped
        and msg_unzipped["jtm_host_name"]
    ):
        db = DbSqlMysql(config=CONFIG)
        # Try to select "pool_name" in workers table by hostname + uselifeLeftrname + poolname
        # Check timediff(now()-end_datetime) in workers table to filter out dead worker
        new_pool_name = jtm_inner_request_q + "." + msg_unzipped["task_pool"]
        num_live_worker_in_pool = db.selectScalar(
            JTM_SQL["select_count_workers_by_jtm_host_name_poolname_enddate"]
            % dict(
                jtm_host_name=msg_unzipped["jtm_host_name"],
                pool_name=new_pool_name,
                hbinterval=w_int * 3,
            ),
            debug=False,
        )
        logger.debug("node cnt in the pool: %s" % num_live_worker_in_pool)

        num_total_num_workers = db.selectScalar(
            JTM_SQL["select_sum_nwpn_workers_by_jtm_host_name_poolname_enddate"]
            % dict(
                jtm_host_name=msg_unzipped["jtm_host_name"],
                pool_name=new_pool_name,
                hbinterval=w_int * 3,
            ),
            debug=False,
        )
        logger.debug("worker cnt in the pool: %s" % num_total_num_workers)

        db.close()
        ret_val = num_total_num_workers if num_total_num_workers else 0
    elif "task_pool" in msg_unzipped and msg_unzipped["task_pool"]:
        db = DbSqlMysql(config=CONFIG)
        # Try to select "pool_name" in workers table by hostname + uselifeLeftrname + poolname
        # Check timediff(now()-end_datetime) in workers table to filter out dead worker
        new_pool_name = jtm_inner_request_q + "." + msg_unzipped["task_pool"]
        num_live_worker_in_pool = db.selectScalar(
            JTM_SQL["select_count_workers_by_poolname_enddate"]
            % dict(pool_name=new_pool_name, hbinterval=w_int * 3),
            debug=False,
        )
        logger.debug("node cnt in the pool: %s" % num_live_worker_in_pool)

        num_total_num_workers = db.selectScalar(
            JTM_SQL["select_sum_nwpn_workers_by_poolname_enddate"]
            % dict(pool_name=new_pool_name, hbinterval=w_int * 3),
            debug=False,
        )
        logger.debug("worker cnt in the pool: %s" % num_total_num_workers)
        db.close()

        ret_val = num_total_num_workers if num_total_num_workers else 0
    elif "jtm_host_name" in msg_unzipped and msg_unzipped["jtm_host_name"]:
        db = DbSqlMysql(config=CONFIG)
        # Check timediff(now() - end_datetime) in workers table to filter out dead worker
        num_live_workers = db.selectScalar(
            JTM_SQL["select_count_workers_by_jtm_host_name"]
            % dict(jtm_host_name=msg_unzipped["jtm_host_name"], hbinterval=w_int * 3),
            debug=False,
        )
        logger.debug("node cnt in the host: %s" % num_live_workers)
        num_total_num_workers = db.selectScalar(
            JTM_SQL["select_sum_nwpn_workers_by_jtm_host_name_enddate"]
            % dict(jtm_host_name=msg_unzipped["jtm_host_name"], hbinterval=w_int * 3),
            debug=False,
        )
        logger.debug("worker cnt in the host: %s" % num_total_num_workers)
        db.close()

        ret_val = num_total_num_workers if num_total_num_workers else 0
    else:
        ret_val = NUM_TOTAL_WORKERS.value

    return ret_val


# -------------------------------------------------------------------------------
def process_remove_pool(task_pool_name):
    """

    :param task_pool_name:
    :return:
    """
    db = DbSqlMysql(config=CONFIG)
    # Get all job id
    zombie_slurm_job_id_list = db.selectAll(
        JTM_SQL["select_all_jid_workers_by_poolname"]
        % dict(
            pool_name=task_pool_name,
        ),
        debug=False,
    )
    db.close()
    logger.debug("zombie_slurm_job_id_list: {}".format(zombie_slurm_job_id_list))

    if len(zombie_slurm_job_id_list) > 0:
        for jid in zombie_slurm_job_id_list:
            scancel_cmd = "scancel %s" % (jid[0])
            _, _, ec = run_sh_command(scancel_cmd, log=logger)
            if ec == 0:
                logger.info("Successfully cancel the job, %s" % (jid[0]))
            else:
                logger.debug("%s not found." % (jid[0]))

    logger.debug("process_remove_pool ****************************************** ")
    db = DbSqlMysql(config=CONFIG)
    db.execute(
        JTM_SQL["update_runs_status_cancelled_by_status_workerId2_poolname"]
        % dict(
            pool_name=task_pool_name,
        ),
        debug=DEBUG,
    )
    # Delete all the workers with pool name from workers table
    db.execute(
        JTM_SQL["delete_from_workers_by_poolname"]
        % dict(pool_name=task_pool_name, pool_name_like=task_pool_name + "_%"),
        debug=DEBUG,
    )
    db.commit()
    db.close()

    logger.info("Pool, %s removed!" % task_pool_name)

    return 1


# -------------------------------------------------------------------------------
def task_kill_proc():
    """
    Check "runs" table if any row with canceled request (canceled=1).
    Check if canceled=1 and wid!=None and status!=-4(which means it's in running status)
    Then, send a kill request to the worker with task ID.

    """
    while True:
        try:
            db = DbSqlMysql(config=CONFIG)
            # Get a list of task ids where canceled -> requested but status != terminated
            tids = [
                int(i)
                for i in db.selectAs1Col(
                    JTM_SQL["select_tids_runs_by_cancelled_and_wid2"]
                )
            ]
            db.close()
        except Exception as e:
            logger.exception(f"task_kill_proc exception in db operation: {e}")

        for tid in tids:
            try:
                # Get the wid and child pid
                db = DbSqlMysql(config=CONFIG)
                worker_id_list = db.selectAs1Col(
                    JTM_SQL["select_workerid_workers_by_tid"] % dict(task_id=tid),
                    debug=False,
                )
                logger.debug("worker list = {}".format(worker_id_list))
                child_proc_id = db.selectScalar(
                    JTM_SQL["select_chilpid_runs_by_tid"] % dict(task_id=tid),
                    debug=False,
                )
                logger.debug("child_proc_id = {}".format(child_proc_id))
                db.close()

                assert child_proc_id > 0
                assert worker_id_list is not None

                # Send task id and process id to worker id
                logger.info(
                    "Send task kill command: {} {} {}".format(
                        tid, worker_id_list, child_proc_id
                    )
                )
                for wid in worker_id_list:
                    send_task_kill_request(tid, wid, child_proc_id)

                # Update runs table with "terminated" status
                logger.debug("Update runs table with terminated for tid {}".format(tid))
                db = DbSqlMysql(config=CONFIG)
                # Here we set cancelled to 2 so that the canceled tids won't be captured
                # again by select_tids_runs_by_cancelled_and_wid2
                db.execute(
                    JTM_SQL["update_runs_status_cancelled_to_terminated_by_tid"]
                    % dict(
                        task_id=tid,
                        status_id=TASK_STATUS["terminated"],
                        now=time.strftime("%Y-%m-%d %H:%M:%S"),
                    ),
                    debug=DEBUG,
                )
                db.commit()
                db.close()

                if not STANDALONE:
                    logger.debug("task kill proc -->")
                    logger.debug(f"status change msg {tid}: '' => terminated")
                    send_update_task_status_msg(
                        tid,
                        None,
                        TASK_STATUS["failed"],
                        fail_code=TASK_STATUS["terminated"],
                    )

            except IndexError as e:
                logger.exception(
                    "Failed to call send_task_kill_request(): {}".format(e)
                )
                logger.debug(
                    "tid, worker_id_list, child_proc_id: {} {} {}".format(
                        tid, worker_id_list, child_proc_id
                    )
                )
            except Exception as e:
                logger.exception(
                    "Something goes wrong in task_kill_proc(): {}".format(e)
                )
                raise

        time.sleep(CONFIG.configparser.getfloat("JTM", "task_kill_interval"))


# -------------------------------------------------------------------------------
def slurm_worker_cleanup_proc():
    """
    Try to find invalid slurm job
    If any, update workers table

    Todo: Very slurm dependent. Do we need this?

    """
    while True:
        try:
            db = DbSqlMysql(config=CONFIG)
            slurm_job_id = [
                int(i)
                for i in db.selectAs1Col(JTM_SQL["select_slurmjid_workers_by_lifeleft"])
            ]
            db.close()
            if len(slurm_job_id) > 0:
                logger.info("Worker checking for %s" % str(slurm_job_id))
            for j in slurm_job_id:
                # Note: slurm dependent code!
                cmd = "squeue -j %d" % j
                so, _, _ = run_sh_command(cmd, log=logger, show_stdout=False)
                if not re.search(r"\bPD\b", so) and not re.search(r"\bR\b", so):
                    logger.info("Found dead worker(s) from %s, %s" % (str(j), so))
                    logger.debug(
                        "slurm_worker_cleanup_proc ****************************** "
                    )
                    db = DbSqlMysql(config=CONFIG)
                    # Note: update status from runs where workerId2 = (select workerId2 from workers)
                    db.execute(
                        JTM_SQL["update_runs_status_cancelled_by_status_workerId2_jid"]
                        % dict(
                            slurm_jid=j,
                        ),
                        debug=DEBUG,
                    )
                    # Delete workers
                    db.execute(
                        JTM_SQL["delete_from_workers_by_slurmjid"]
                        % dict(
                            slurm_jid=j,
                        ),
                        debug=DEBUG,
                    )
                    db.commit()
                    db.close()
                else:
                    logger.debug("Job id {} is alive.".format(j))

                time.sleep(1)

            time.sleep(CONFIG.configparser.getfloat("JTM", "worker_kill_interval"))
        except Exception as e:
            logger.exception(f"slurm_worker_cleanup_proc exception: {e}")
            raise


# -------------------------------------------------------------------------------
def proc_clean_exit(plist):
    """
    Kill all parent + child processes and clean them all
    WHEN one of child processes is abnormally crashed.
    And, when new deployment, supervisord will send SIGTERM to
    the JTM manager.

    :param plist: process handle list
    :return:
    """
    for p in plist:
        try:
            if p is not None and p.is_alive():
                p.terminate()
        except AssertionError:
            # is_alive() raises AssertionError
            # if assert self._parent_pid != os.getpid()
            logger.warning("Skipping is_alive() checking for the parent process.")
        except Exception as e:
            # print log and just pass
            logger.exception(f"Failed to terminate a child process: {e}")

    os._exit(1)


# -------------------------------------------------------------------------------
def kill_child_proc(ppid):
    try:
        for process in psutil.process_iter():
            _ppid = process.ppid()
            if _ppid == ppid:
                _pid = process.pid
                if sys.platform == "win32":
                    process.terminate()
                else:
                    os.system("kill -9 {0}".format(_pid))
    except Exception as e:
        logger.exception(f"Failed to terminate a child process: {e}")
        raise


# -------------------------------------------------------------------------------
def check_num_threads(mode: str, n_manager_threads: int) -> bool:
    ps_cmd = "ps -aef | grep -v grep | grep -v check-manager | grep jtm | grep manager "
    if mode != "test":
        ps_cmd += f"| grep jaws-{mode} "
    else:
        ps_cmd += "| grep test "
    ps_cmd += "| wc -l"
    num_total_procs = 0
    try:
        so, _, ec = run_sh_command(ps_cmd, log=logger, show_stdout=False)
        num_total_procs = int(so.rstrip())
    except TypeError as te:
        logger.exception(te)
        logger.error(so)
        return False
    except Exception as e:
        logger.exception(e)
        logger.error(ps_cmd)
        return False
    else:
        # THIS IS THE NUM OF JTM MANAGER PROCS = 8
        if num_total_procs != n_manager_threads:
            return False
        else:
            return True


# -------------------------------------------------------------------------------
def manager(
    ctx: object, custom_log_dir_name: str, b_resource_usage_log_on: bool
) -> int:
    """

    :param ctx:
    :param custom_log_dir_name:
    :param b_resource_usage_log_on:
    :return:
    """
    global CONFIG
    CONFIG = ctx.obj["config"]
    global DEBUG
    DEBUG = ctx.obj["debug"]
    # config file has precedence
    config_debug = CONFIG.configparser.getboolean("SITE", "debug")
    if config_debug:
        DEBUG = config_debug
    global TASK_STATUS
    TASK_STATUS = CONFIG.constants.TASK_STATUS
    jtm_task_request_q = CONFIG.configparser.get("JTM", "jtm_task_request_q")

    # if STANDALONE == 1, send_update_task_status_msg(), which sends task status change
    # RMQ message to JAWS Site, won't be called
    global STANDALONE
    try:
        STANDALONE = CONFIG.configparser.getboolean("JTM", "standalone")
    except Exception:
        STANDALONE = False

    # Log dir setting
    datetime_str = datetime.datetime.now().strftime("%Y-%m-%d")
    log_dir_name = CONFIG.configparser.get("JTM", "log_dir")
    if custom_log_dir_name:
        log_dir_name = custom_log_dir_name
    log_dir_name = os.path.join(log_dir_name, "manager")
    make_dir(log_dir_name)
    log_file_name = "%s/jtm_%s.log" % (log_dir_name, datetime_str)

    log_level = "info"
    if DEBUG:
        log_level = "debug"

    setup_custom_logger(log_level, log_dir_name, log_file_name, 1, 1)

    logger.info("JTM Manager, version: {}".format(CONFIG.constants.VERSION))
    logger.info(
        "\n*****************\nDebug mode is %s\n*****************"
        % ("ON" if DEBUG else "OFF")
    )
    logger.info("Set jtm log file location to %s", log_dir_name)

    prod_mod = False
    if CONFIG.configparser.get("JTM", "run_mode") == "prod":
        prod_mod = True

    logger.info(
        "JTM main exchange: %s", CONFIG.configparser.get("JTM", "jgi_jtm_main_exch")
    )
    logger.info("JTM config file: %s", CONFIG.config_file)
    logger.info("RabbitMQ broker: %s", CONFIG.configparser.get("RMQ", "host"))
    logger.info("Default task queue name: %s", jtm_task_request_q)
    logger.info(
        "Default result queue name: %s",
        CONFIG.configparser.get("JTM", "jtm_task_result_q"),
    )
    logger.info("Database server: %s", CONFIG.configparser.get("MYSQL", "host"))
    logger.info("Database user name: %s", CONFIG.configparser.get("MYSQL", "user"))
    logger.info("Database name: %s", CONFIG.configparser.get("MYSQL", "db"))
    logger.info("JTM user name: %s", CONFIG.configparser.get("SITE", "user_name"))
    logger.info(
        "\n*****************\nRun mode is %s\n*****************"
        % ("PROD" if prod_mod else "DEV")
    )

    # MySQL: prepare task table
    db = DbSqlMysql(config=CONFIG)
    db.ddl(JTM_SQL["create_database"] % CONFIG.configparser.get("MYSQL", "db"))
    db.ddl(JTM_SQL["use_database"] % CONFIG.configparser.get("MYSQL", "db"))
    db.ddl(JTM_SQL["create_table_tasks"])
    db.ddl(JTM_SQL["create_table_runs"])
    db.ddl(JTM_SQL["create_table_workers"])
    db.close()

    # pid list
    plist = list()

    def signal_handler(signum, frame):
        proc_clean_exit(plist)

    signal.signal(signal.SIGTERM, signal_handler)
    ppid = os.getpid()
    logger.debug(f"main_proc pid = {ppid}")

    # Start heartbeat receiving proc
    worker_hb_queue_name = CONFIG.configparser.get("JTM", "worker_hb_q_postfix")
    try:
        recv_hb_from_worker_proc_hdl = mp.Process(
            target=recv_hb_from_worker_proc,
            args=(worker_hb_queue_name, log_dir_name, b_resource_usage_log_on),
        )
        recv_hb_from_worker_proc_hdl.start()
        plist.append(recv_hb_from_worker_proc_hdl)
        logger.debug(
            f"recv_hb_from_worker_proc pid = {recv_hb_from_worker_proc_hdl.pid}"
        )
    except Exception as e:
        logger.exception("recv_hb_from_worker_proc: {}".format(e))
        proc_clean_exit(plist)
        raise

    # Start cancelled task cleanup proc
    try:
        task_kill_proc_hdl = mp.Process(target=task_kill_proc)
        task_kill_proc_hdl.start()
        plist.append(task_kill_proc_hdl)
        logger.debug(f"task_kill_proc pid = {task_kill_proc_hdl.pid}")
    except Exception as e:
        logger.exception("task_kill_proc: {}".format(e))
        proc_clean_exit(plist)
        raise

    # Start worker cleanup proc
    try:
        worker_cleanup_proc_hdl = mp.Process(target=slurm_worker_cleanup_proc)
        worker_cleanup_proc_hdl.start()
        plist.append(worker_cleanup_proc_hdl)
        logger.debug(f"worker_cleanup_proc pid = {worker_cleanup_proc_hdl.pid}")
    except Exception as e:
        logger.exception("slurm_worker_cleanup_proc: {}".format(e))
        proc_clean_exit(plist)
        raise

    # Start result receiving proc
    try:
        recv_result_from_worker_proc_hdl = mp.Process(
            target=WorkerResultReceiver(config=CONFIG).start
        )
        recv_result_from_worker_proc_hdl.start()
        plist.append(recv_result_from_worker_proc_hdl)
        logger.debug(
            f"recv_result_from_worker_proc pid = {recv_result_from_worker_proc_hdl.pid}"
        )
    except Exception as e:
        logger.exception("WorkerResultReceiver: {}".format(e))
        proc_clean_exit(plist)
        raise

    # Start JtmCommandRunner
    try:
        process_jtminterface_command_proc_hdl = mp.Process(
            target=JtmCommandRunner(config=CONFIG, ppid=ppid).start
        )
        process_jtminterface_command_proc_hdl.start()
        plist.append(process_jtminterface_command_proc_hdl)
        logger.debug(
            f"process_jtminterface_taskcommand_proc pid = {process_jtminterface_command_proc_hdl.pid}"
        )
    except Exception as e:
        logger.exception("JtmCommandRunner: {}".format(e))
        proc_clean_exit(plist)
        raise

    # Start JtmNonTaskCommandRunner
    try:
        process_jtminterface_nontaskcommand_proc_hdl = mp.Process(
            target=JtmNonTaskCommandRunner(config=CONFIG, ppid=ppid).start
        )
        process_jtminterface_nontaskcommand_proc_hdl.start()
        plist.append(process_jtminterface_nontaskcommand_proc_hdl)
        logger.debug(
            f"process_jtminterface_nontaskcommand_proc pid = {process_jtminterface_nontaskcommand_proc_hdl.pid}"
        )
    except Exception as e:
        logger.exception("JtmNonTaskCommandRunner: {}".format(e))
        proc_clean_exit(plist)
        raise

    # Start JTM JSON-RPC server for monitoring
    def jtm_manager_status(params):
        # 'params' is not used
        alive = True
        run_mode = CONFIG.configparser.get("JTM", "run_mode")
        n_manager_threads = CONFIG.constants.JTM_NUM_PROCS
        if not check_num_threads(run_mode, n_manager_threads):
            alive = False
        if alive:
            return responses.success(True)
        else:
            return responses.success(False)

    operations = {
        "server_status": {
            "function": jtm_manager_status,
            "required_params": [],
        }
    }
    jtm_rpc_server_params = CONFIG.configparser._sections("JTM_RPC_SERVER")
    logger.debug("jtm_rpc_server params: %s", jtm_rpc_server_params)
    app = rpc_server.RpcServer(jtm_rpc_server_params, operations)
    app.start_server()

    logger.info("Waiting for worker's heartbeats from %s", worker_hb_queue_name)
    logger.info("Waiting for a task request from %s", jtm_task_request_q)


    return 0
