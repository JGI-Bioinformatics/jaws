#! /usr/bin/env python
# pylint: disable=C0111,C0103,R0205
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)

"""

jtm worker


Example scenario
1. jtm_submit sends a msg to "jgi_main_exchange" with "jtm_task_request_queue" tag.
2. jtm listens to "jtm_task_request_queue" which is bound to "jgi_main_exchange"
3. when a task is ready, jtm takes it and sends it to one of the workers
   (to jgi_jtm_inner_main_exchange)
4. workers listen to jgi_jtm_inner_request_queue which is bound to "jgi_jtm_inner_main_exchange"
5. when a task is done, a worker sends a result msg to "jgi_jtm_inner_main_exchange" with
   "jgi_jtm_inner_result_queue" tag
6. jtm listens to "jgi_jtm_inner_result_queue" queue; when a result is ready, takes and updates
   tables

"""
import datetime
import json
import multiprocessing as mp
import os
import resource
import signal
import socket
import subprocess
import sys
import time
import amqpstorm
import psutil
import shortuuid

from jaws_jtm.common import logger, setup_custom_logger
from jaws_jtm.lib.msgcompress import zdumps, zloads
from jaws_jtm.lib.rabbitmqconnection import JtmAmqpstormBase, RmqConnectionAmqpstorm
from jaws_jtm.lib.resourceusage import (
    get_cpu_load,
    get_free_memory,
    get_pid_tree,
    get_resident_memory_usage,
    get_runtime,
    get_total_mem_usage_per_node,
    get_virtual_memory_usage,
)
from jaws_jtm.lib.run import make_dir, run_sh_command

# -------------------------------------------------------------------------------
# Globals
# -------------------------------------------------------------------------------
# This ipc pipe is to send task_id (by do_work())to send_hb_to_client_proc()
# when a task is requested
PIPE_TASK_ID_SEND, PIPE_TASK_ID_RECV = mp.Pipe()
# To share user task proc id
USER_PROC_PROC_ID = mp.Value("i", 0)
# To share this worker life left
WORKER_LIFE_LEFT_IN_MINUTE = mp.Value("i", 0)
UNIQ_WORKER_ID = str(shortuuid.uuid())
WORKER_START_TIME = datetime.datetime.now()
PARENT_PROCESS_ID = os.getpid()  # parent process id
THIS_WORKER_TYPE = None


class TaskTerminator(JtmAmqpstormBase):
    """
    Process task cancellation

    """

    def start(self):
        """Start the OnWorkerResult.
        :return:
        """
        if not self.connection:
            self.create_connection()
        while True:
            try:
                channel = self.connection.channel()
                channel.basic.qos(prefetch_count=1)
                channel.exchange.declare(
                    exchange=self.jtm_task_kill_exch,
                    exchange_type="fanout",
                    durable=True,
                    auto_delete=False,
                )
                channel.queue.declare(
                    queue=self.jtm_task_kill_q,
                    durable=True,
                    exclusive=False,
                    auto_delete=True,
                )
                channel.queue.bind(
                    exchange=self.jtm_task_kill_exch,
                    queue=self.jtm_task_kill_q,
                    routing_key=UNIQ_WORKER_ID,
                )
                channel.basic.consume(
                    self.process_kill, self.jtm_task_kill_q, no_ack=False
                )
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
                raise

    def process_kill(self, message):
        """
        Callback for processing user task

        :param message:
        :return:
        """
        msg_unzipped = json.loads(zloads(message.body))
        logger.debug("Received task termination command: %r" % msg_unzipped)

        if msg_unzipped["worker_id"] == UNIQ_WORKER_ID:
            if msg_unzipped["child_pid"] > 0:
                logger.info("Process termination request received.")

                # This -9 is to notify run_user_task() that it's killed by user requests
                # Also send_hb_to_client() will check this for adjust childpid and parentpid
                # Note: this should done first to signal run_user_task() that the process is killed.
                USER_PROC_PROC_ID.value = -9

                # kill if there is child's children
                for i in get_pid_tree(msg_unzipped["child_pid"]):
                    kill_cmd = "kill -9 %d" % i
                    logger.info(
                        f"Executing {kill_cmd} for taskID, {msg_unzipped['task_id']}"
                    )
                    _, _, ec = run_sh_command(kill_cmd, log=logger)
                    if ec == 0:
                        logger.info("Successfully terminate a user task process.")
                    else:
                        logger.warning(
                            f"User process not found. Ignore the termination command, {kill_cmd}"
                        )

                # Kill the main child process
                # Note: can consider to use "pkill -9 -P ppid" to kill the family
                kill_cmd = f"kill -9 {msg_unzipped['child_pid']}"
                logger.info(
                    f"Executing {kill_cmd} for taskID, {msg_unzipped['task_id']}"
                )
                _, _, ec = run_sh_command(kill_cmd, log=logger)
                if ec == 0:
                    logger.info("Successfully terminate a user task process.")
                else:
                    logger.warning(
                        f"User process not found. Failed to execute the command, {kill_cmd}"
                    )
            else:
                logger.warning("No valid child process id to terminate.")

            message.ack()

        else:
            # If UNIQ_WORKER_ID is not for me, reject and requeue it
            logger.info("Send NACK to the broker")
            message.nack(requeue=True)


class TaskRunner(JtmAmqpstormBase):
    """
    Process user task

    """

    def start(self, req_q):
        """Start the OnWorkerResult.
        :return:
        """
        if not self.connection:
            self.create_connection()
        while True:
            try:
                channel = self.connection.channel()
                channel.basic.qos(prefetch_count=1)
                channel.exchange.declare(
                    exchange=self.jtm_inner_main_exch,
                    exchange_type="direct",
                    durable=True,
                    auto_delete=False,
                )
                channel.queue.declare(
                    queue=req_q, durable=True, exclusive=False, auto_delete=True
                )
                channel.queue.bind(
                    exchange=self.jtm_inner_main_exch, queue=req_q, routing_key=req_q
                )
                channel.basic.consume(self.process_task, req_q, no_ack=False)
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
                raise

    def process_task(self, message):
        """
        Callback for processing user task

        :param message:
        :return:
        """
        msg_unzipped = json.loads(zloads(message.body))
        task_id = msg_unzipped["task_id"]

        # Send taskid to send_hb_to_client_proc() so that the task_id can be shown in the hb messages
        PIPE_TASK_ID_SEND.send(task_id)

        logger.info("Received a task, %r" % (msg_unzipped,))
        logger.debug(f"Return queue = {message.reply_to}")

        result_dict = run_user_task(msg_unzipped)

        json_data = json.dumps(result_dict)
        logger.debug(f"Reply msg with result: {str(json_data)}")
        logger.debug(json_data)
        response = zdumps(json_data)
        logger.debug("Send ACK to the manager")

        logger.debug("Send result to the manager")
        try:
            with RmqConnectionAmqpstorm(config=CONFIG).open() as conn:
                with conn.channel() as ch:
                    ch.exchange.declare(
                        exchange=self.jtm_inner_main_exch,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False,
                    )
                    ch.queue.declare(
                        queue=self.inner_result_queue_name,
                        durable=True,
                        exclusive=False,
                        auto_delete=True,
                    )
                    ch.queue.bind(
                        exchange=self.jtm_inner_main_exch,
                        queue=self.inner_result_queue_name,
                        routing_key=self.inner_result_queue_name,
                    )
                    properties = {
                        "correlation_id": message.correlation_id,
                        "reply_to": message.reply_to,
                    }
                    response = amqpstorm.Message.create(ch, response, properties)
                    response.publish(self.inner_result_queue_name)
        except Exception as detail:
            logger.exception(
                f"Exception: Failed to send a request to {self.inner_result_queue_name}"
            )
            logger.exception(f"Detail: {detail}")
            raise OSError(2, "Failed to send a request to a worker")

        try:
            message.ack()
        except Exception as detail:
            logger.exception(
                "Exception: Failed to send an ACK for the task msg"
            )
            logger.exception(f"Detail: {detail}")
            PIPE_TASK_ID_SEND.send(0)
            USER_PROC_PROC_ID.value = 0
            raise OSError(2, "Failed to send an ACK for the task completed")

        # Send taskid=0 to send_hb_to_client_thread() b/c the requested task is completed
        PIPE_TASK_ID_SEND.send(0)

        # Reset child pid
        USER_PROC_PROC_ID.value = 0


# -------------------------------------------------------------------------------
def run_user_task(msg_unzipped):
    """
    Run a user command in msg_zipped_to_send

    :param msg_unzipped: uncompressed msg from client
    :return:
    """
    return_msg = {}

    # Uncompress msg to get a task
    logger.info(f"User task to run: {msg_unzipped}")
    task_id = msg_unzipped["task_id"]
    user_task_cmd = msg_unzipped["user_cmd"]
    task_type = msg_unzipped["task_type"]

    return_msg["task_id"] = task_id
    return_msg["user_cmd"] = user_task_cmd
    return_msg["task_type"] = task_type

    if "stdout" in msg_unzipped:
        user_task_cmd + f" > {msg_unzipped['stdout']}"
    if "stderr" in msg_unzipped:
        user_task_cmd + f" 2>{msg_unzipped['stderr']}"

    # Run the task
    logger.info(f"Running task ID {task_id}...")

    p = None
    time_out_in_minute = 0
    done_f = CONFIG.constants.DONE_FLAGS
    w_int = CONFIG.configparser.getfloat("JTM", "worker_hb_recv_interval")

    if THIS_WORKER_TYPE != "manual":
        # Wait until WORKER_LIFE_LEFT_IN_MINUTE is updated by worker's HB
        limit = 0
        if WORKER_LIFE_LEFT_IN_MINUTE.value <= 0:
            while True:
                logger.debug(f"worker life: {WORKER_LIFE_LEFT_IN_MINUTE.value}")
                if WORKER_LIFE_LEFT_IN_MINUTE.value > 0:
                    break
                if limit == w_int * 3:
                    logger.info("Worker timeout.")
                    return_msg["done_flag"] = str(done_f["failed with timeout"])
                    return_msg["ret_msg"] = "User task timeout"
                    return return_msg
                time.sleep(1)
                limit += 1

        # ex) WORKER_LIFE_LEFT_IN_MINUTE = 20min and TASK_KILL_TIMEOUT_MINUTE = 10min
        # timeout will be set as 10min
        # TASK_KILL_TIMEOUT_MINUTE is a extra housekeeping time after explicitly
        # terminate a task.
        tkill_time = CONFIG.configparser.getint("JTM", "task_kill_timeout_minute")
        assert tkill_time
        logger.debug(f"worker life: {WORKER_LIFE_LEFT_IN_MINUTE.value}")
        time_out_in_minute = int(WORKER_LIFE_LEFT_IN_MINUTE.value - tkill_time)
        logger.info(f"Timeout in minute: {time_out_in_minute}")

        if time_out_in_minute <= tkill_time:
            logger.info(f"Not enough time to run the task, {task_id}.")
            return_msg["done_flag"] = str(done_f["failed with timeout"])
            return_msg["ret_msg"] = "User task timeout"
            return return_msg

    proc_return_code = -1

    logger.debug("Start subprocess to run a task.")
    try:
        p = subprocess.Popen(
            user_task_cmd.split(),
            env=os.environ,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    except MemoryError:
        logger.exception(f"Exception: Out of memory, {user_task_cmd}")
        proc_return_code = done_f["failed with out-of-mem"]
    except FileNotFoundError:
        logger.exception(f"Exception: Input file or command not found, {user_task_cmd}")
        proc_return_code = 2
    except Exception as detail:
        logger.exception(
            f"Exception: Failed to run user command, {user_task_cmd}",
        )
        logger.exception(f"Detail: {detail}")
    else:
        # Set USER_PROC_PROC_ID = forked child process id (this value will be sent
        # to send_hb_to_client function)
        USER_PROC_PROC_ID.value = p.pid
        logger.debug(f"USER_PROC_PROC_ID.value {USER_PROC_PROC_ID.value}")
        logger.debug(f"User command: {user_task_cmd}")

    stdout_str = None
    time_out_in_second = None
    if time_out_in_minute != 0:
        time_out_in_second = time_out_in_minute * 60

    if p is not None:
        # ref)
        # https://www.programcreek.com/python/example/56781/subprocess.TimeoutExpired
        try:
            stdout_str = p.communicate(timeout=time_out_in_second)[0]
            if type(stdout_str) == bytes:
                stdout_str = stdout_str.decode()
        except subprocess.TimeoutExpired:
            logger.exception("subprocess call timeout")
            p.kill()
            proc_return_code = done_f["failed with timeout"]
        else:
            p.wait()
            p.poll()
            if stdout_str.find("command not found") != -1:
                proc_return_code = 2  # set it as command-not-found error
            else:
                proc_return_code = p.returncode
    else:
        logger.error("subprocess call failed")

    # Prepare result to send back
    return_msg["worker_id"] = UNIQ_WORKER_ID
    return_msg["host_name"] = socket.gethostname()

    if USER_PROC_PROC_ID.value == -9:
        return_msg["done_flag"] = done_f["failed with user termination"]
        return_msg["ret_msg"] = "Task cancelled."
    else:
        if proc_return_code == done_f["failed with timeout"]:
            logger.info("User task timeout. Not enough time left for the worker.")
            return_msg["done_flag"] = str(done_f["failed with timeout"])
            return_msg["ret_msg"] = "User task timeout"
        elif proc_return_code == done_f["failed with out-of-mem"]:
            logger.info("User task out-of-mem.")
            return_msg["done_flag"] = str(done_f["failed with out-of-mem"])
            return_msg["ret_msg"] = "User task out-of-mem"
        elif proc_return_code == 2:  # system code
            logger.info("input file or command not found.")
            return_msg["done_flag"] = done_f[
                "failed with input file or command not found"
            ]
            return_msg["ret_msg"] = "Input file or command not found."
        elif proc_return_code == 0:  # system code
            logger.info(f"Task# {task_id} completed!")
            return_msg["done_flag"] = done_f["success"]
            return_msg["ret_msg"] = ""
        else:
            logger.critical(
                f"Failed to execute a task, {user_task_cmd}. Non-zero exit code. stdout = {stdout_str}."
            )
            return_msg["done_flag"] = done_f["failed to run user command"]
            return_msg["ret_msg"] = stdout_str

    logger.info(f"Reply msg prepared with result: {return_msg}")

    return return_msg


# -------------------------------------------------------------------------------
def send_hb_to_client_proc(
    interval,
    slurm_job_id,
    mem_per_node,
    mem_per_core,
    num_cores,
    job_time,
    pool_name,
    nwpn,
    show_resource_log,
):
    """
    Send heartbeats to the client

    :param interval: time interval to send heartbeats to the client
    :param slurm_job_id: SLURM job id
    :param mem_per_node: memory request per node
    :param mem_per_core: memory request per core
    :param num_cores: number of cores
    :param job_time: wallclocktime
    :param pool_name: pool name
    :param nwpn: number of workers per node
    :return:
    """
    host_name = socket.gethostname()
    exch_name = CONFIG.configparser.get("JTM", "jtm_worker_hb_exch")
    worker_hb_queue = CONFIG.configparser.get("JTM", "worker_hb_q_postfix")
    today = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    try:
        ip_address = socket.gethostbyname(host_name)
    except Exception:
        ip_address = None
    ncore_param = mp.cpu_count()
    if CONFIG.constants.WORKER_TYPE[THIS_WORKER_TYPE] > 0:
        ncore_param = num_cores
    hb_msg = CONFIG.constants.HB_MSG
    w_type = CONFIG.constants.WORKER_TYPE
    task_id = 0
    proc_id_list_merged = list()
    vmem_usage_list = list()
    rmem_usage_list = list()
    cpu_load_list = list()
    end_date_time = None
    while True:
        # NOTE: must cope with the fast mem consumption at the begining of the process
        if USER_PROC_PROC_ID.value == 0:
            child_pid = PARENT_PROCESS_ID
        else:
            child_pid = int(USER_PROC_PROC_ID.value)

        # Collect pids from process tree
        try:
            root_pid_num = get_pid_tree(PARENT_PROCESS_ID)
            if PARENT_PROCESS_ID != child_pid:
                if child_pid == -9:  # if it's terminated by "kill"
                    child_pid = PARENT_PROCESS_ID
                pid_list_child = get_pid_tree(child_pid)
                if len(pid_list_child) > 0:
                    proc_id_list_merged = root_pid_num + pid_list_child[1:]
                else:
                    # NOTE: be careful on this resetting! Might lose child pid
                    # USER_PROC_PROC_ID.value = 0
                    proc_id_list_merged = root_pid_num
            else:
                proc_id_list_merged = root_pid_num
        except Exception as e:
            logger.warning(f"get_pid_tree error: {e}")

        # Collect vmem usage from process tree
        for pid in proc_id_list_merged:
            try:
                vmem = get_virtual_memory_usage(pid, 0.0, False)
            except ValueError:
                logger.warning("ValueError: Failed to collect VM memory usage.")
            except UnboundLocalError:
                logger.warning("UnboundLocalError: No entry in process id list.")
            else:
                vmem_usage_list.append(vmem)

        # Collect rss mem usage from process tree
        for pid in proc_id_list_merged:
            try:
                rmem = get_resident_memory_usage(pid, 0.0, False)
            except ValueError:
                logger.warning("ValueError: Failed to collect RES memory usage.")
            except UnboundLocalError:
                logger.warning("UnboundLocalError: No entry in process id list.")
            else:
                rmem_usage_list.append(rmem)

        # Collect mem_usages for all pids in the tree and get sum()
        rmem_usage = "%.1f" % sum(rmem_usage_list) if len(rmem_usage_list) > 0 else 0.0
        vmem_usage = "%.1f" % sum(vmem_usage_list) if len(vmem_usage_list) > 0 else 0.0

        # Collect cpu_usages for all pids in the tree and get max()
        for pid in proc_id_list_merged:
            try:
                cload = get_cpu_load(pid)
            except Exception as e:
                logger.warning(f"get_cpu_load() exception: {e}")
            else:
                cpu_load_list.append(cload)

        max_cpu_load = max(cpu_load_list) if len(cpu_load_list) > 0 else 0.0

        # Only get the run time of child_pid
        if sys.platform.lower() == "darwin":
            # Todo: Add a method to get etime on Mac OS
            proc_run_time = ""
        else:
            proc_run_time = get_runtime(child_pid)

        if max_cpu_load == "":
            max_cpu_load = 0.0

        if proc_run_time == "":
            proc_run_time = 0

        # Get % mem used per node
        # This is for node-based scheduling
        try:
            perc_used_mem = "%.1f" % get_total_mem_usage_per_node()
        except Exception as e:
            logger.warning(f"get_total_mem_usage_per_node() exception: {e}")
            perc_used_mem = 0.0

        # Check if there is any task id in the ipc pipe
        if PIPE_TASK_ID_RECV.poll():
            task_id = PIPE_TASK_ID_RECV.recv()

        if slurm_job_id:
            temp = WORKER_LIFE_LEFT_IN_MINUTE.value
            try:
                # hh:mm:ss --> seconds
                job_runtime_in_sec = (
                    int(job_time.split(":")[0]) * 3600
                    + int(job_time.split(":")[1]) * 60
                    + int(job_time.split(":")[2])
                )
                end_date_time = WORKER_START_TIME + datetime.timedelta(
                    seconds=int(job_runtime_in_sec)
                )
                delta = end_date_time - datetime.datetime.now()
                WORKER_LIFE_LEFT_IN_MINUTE.value = int(
                    divmod(delta.total_seconds(), 60)[0]
                )
            except Exception as e:
                logger.warning(
                    f"Something wrong in computing remaining wall clock time: {e}"
                )
                WORKER_LIFE_LEFT_IN_MINUTE.value = temp

        msg_dict_to_send = {
            hb_msg["child_pid"]: child_pid,
            hb_msg["clone_time_rate"]: 0.0,  # OBSOLETE
            hb_msg["cpu_load"]: max_cpu_load,
            hb_msg["end_date"]: today,
            hb_msg["host_name"]: host_name,
            hb_msg["ip_address"]: ip_address,
            hb_msg["job_time"]: job_time if w_type[THIS_WORKER_TYPE] > 0 else None,
            hb_msg["jtm_host_name"]: CONFIG.configparser.get("SITE", "jtm_host_name"),
            hb_msg["life_left"]: WORKER_LIFE_LEFT_IN_MINUTE.value,
            hb_msg["mem_per_core"]: mem_per_core
            if w_type[THIS_WORKER_TYPE] > 0
            else "",
            hb_msg["mem_per_node"]: mem_per_node
            if w_type[THIS_WORKER_TYPE] > 0
            else "",
            hb_msg["num_cores"]: ncore_param,
            hb_msg["num_tasks"]: 0,
            hb_msg["num_workers_on_node"]: 1,  # not used
            hb_msg["perc_mem_used"]: perc_used_mem,
            hb_msg["pool_name"]: pool_name,
            hb_msg["ret_msg"]: "hb",
            hb_msg["rmem_usage"]: rmem_usage,
            hb_msg["root_pid"]: PARENT_PROCESS_ID,
            hb_msg["run_time"]: proc_run_time,
            hb_msg["slurm_jobid"]: slurm_job_id,
            hb_msg["task_id"]: task_id,
            hb_msg["vmem_usage"]: vmem_usage,
            hb_msg["worker_id"]: UNIQ_WORKER_ID,
            hb_msg["worker_type"]: w_type[THIS_WORKER_TYPE],
            hb_msg["nwpn"]: nwpn,  # not used
        }

        msg_zipped_to_send = zdumps(json.dumps(msg_dict_to_send))

        if show_resource_log:
            logger.info(msg_dict_to_send)

        max_retries = CONFIG.configparser.getint("JTM", "max_retries")
        exp = CONFIG.configparser.get("JTM", "worker_s_hb_expiration")
        attempts = 0
        while True:
            attempts += 1
            with RmqConnectionAmqpstorm(config=CONFIG).open() as conn:
                with conn.channel() as ch:
                    ch.exchange.declare(
                        exchange=exch_name,
                        exchange_type="direct",
                        durable=False,
                        auto_delete=False,
                    )
                    message = amqpstorm.Message.create(
                        ch, msg_zipped_to_send, properties={"expiration": exp}
                    )
                    try:
                        message.publish(worker_hb_queue)
                        break
                    except amqpstorm.AMQPError as why:
                        logger.exception(why)
                        logger.warning(
                            f"Retry to send a heartbeat to the manager: trial #{attempts}"
                        )
                        if attempts > max_retries:
                            logger.critical(
                                f"Failed to send a heartbeat to the manageer: message = {msg_dict_to_send}"
                            )
                            raise amqpstorm.AMQPError
                        time.sleep(min(attempts * 2, 30))
                    except Exception as e:
                        logger.exception(e)
                        raise

        time.sleep(interval)


# -------------------------------------------------------------------------------
def proc_clean_exit(pid_list):
    """

    :param pid_list: process handle list
    :return:
    """
    for p in pid_list:
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
def worker(
    ctx: object,
    heartbeat_interval_param: int,
    custom_log_dir: str,
    custom_job_log_dir_name: str,
    pool_name_param: str,
    dry_run: bool,
    slurm_job_id_param: int,
    worker_type_param: str,
    cluster_name_param: str,
    worker_clone_time_rate_param: float,
    num_workers_per_node_param: int,
    worker_id_param: str,
    charging_account_param: str,
    num_nodes_to_request_param: int,
    num_cores_to_request_param: int,
    constraint_param: str,
    mem_per_node_to_request_param: str,
    mem_per_cpu_to_request_param: str,
    qos_param: str,
    partition_param: str,
    job_time_to_request_param: str,
    show_resource_log: bool,
) -> int:
    """

    :param ctx:
    :param heartbeat_interval_param:
    :param custom_log_dir:
    :param custom_job_log_dir_name:
    :param pool_name_param:
    :param dry_run:
    :param slurm_job_id_param:
    :param worker_type_param:
    :param cluster_name_param:
    :param worker_clone_time_rate_param:
    :param num_workers_per_node_param:
    :param worker_id_param:
    :param charging_account_param:
    :param num_nodes_to_request_param:
    :param num_cores_to_request_param:
    :param constraint_param:
    :param mem_per_node_to_request_param:
    :param mem_per_cpu_to_request_param:
    :param qos_param:
    :param partition_param:
    :param job_time_to_request_param:
    :param show_resource_log:
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
    global JTM_INNER_MAIN_EXCH
    JTM_INNER_MAIN_EXCH = CONFIG.configparser.get("JTM", "jtm_inner_main_exch")

    prod_mod = False
    if CONFIG.configparser.get("JTM", "run_mode") == "prod":
        prod_mod = True

    # Set uniq worker id if worker id is provided in the params
    if worker_id_param:
        global UNIQ_WORKER_ID
        UNIQ_WORKER_ID = worker_id_param

    # Job log dir setting
    job_script_dir_name = os.path.join(CONFIG.configparser.get("JTM", "log_dir"), "job")
    if custom_job_log_dir_name:
        job_script_dir_name = custom_job_log_dir_name
    make_dir(job_script_dir_name)

    # Log dir setting
    datetime_str = datetime.datetime.now().strftime("%Y-%m-%d")
    datetime_str = UNIQ_WORKER_ID + "_" + datetime_str
    log_dir_name = CONFIG.configparser.get("JTM", "log_dir")
    if custom_log_dir:
        log_dir_name = custom_log_dir
    log_dir_name = os.path.join(log_dir_name, "worker")
    make_dir(log_dir_name)
    log_file_name = f"{log_dir_name}/jtm_worker_{datetime_str}.log"

    # Logger setting
    log_level = "info"
    if DEBUG:
        log_level = "debug"

    setup_custom_logger(log_level, log_dir_name, log_file_name, 1, 1)

    logger.info(f"JTM Worker, version: {CONFIG.constants.VERSION}")
    hearbeat_interval = CONFIG.configparser.getfloat("JTM", "worker_hb_send_interval")
    logger.info(
        "\n*****************\nDebug mode is %s\n*****************"
        % ("ON" if DEBUG else "OFF")
    )
    logger.info("Set jtm log file location to %s", log_dir_name)
    logger.info("Set jtm job file location to %s", job_script_dir_name)
    logger.info("RabbitMQ broker: %s", CONFIG.configparser.get("RMQ", "host"))
    logger.info("RabbitMQ port: %s", CONFIG.configparser.get("RMQ", "port"))
    logger.info("JTM user name: %s", CONFIG.configparser.get("SITE", "user_name"))
    logger.info("Unique worker ID: %s", UNIQ_WORKER_ID)
    logger.info(
        "\n*****************\nRun mode is %s\n*****************"
        % ("PROD" if prod_mod else "DEV")
    )
    logger.info("env activation: %s", CONFIG.configparser.get("JTM", "env_activation"))
    logger.info("JTM config file: %s" % (CONFIG.config_file))

    # Slurm config
    num_nodes_to_request = (
        num_nodes_to_request_param
        if num_nodes_to_request_param
        else CONFIG.configparser.getint("SLURM", "nnodes")
    )
    num_workers_per_node = (
        num_workers_per_node_param
        if num_workers_per_node_param
        else CONFIG.configparser.getint("JTM", "num_workers_per_node")
    )
    assert num_workers_per_node > 0
    mem_per_cpu_to_request = (
        mem_per_cpu_to_request_param
        if mem_per_cpu_to_request_param
        else CONFIG.configparser.get("SLURM", "mempercpu")
    )
    mem_per_node_to_request = (
        mem_per_node_to_request_param
        if mem_per_node_to_request_param
        else CONFIG.configparser.get("SLURM", "mempernode")
    )
    assert mem_per_cpu_to_request
    assert mem_per_node_to_request
    num_cpus_to_request = (
        num_cores_to_request_param
        if num_cores_to_request_param
        else CONFIG.configparser.getint("SLURM", "ncpus")
    )
    assert num_cpus_to_request

    # Set CPU affinity for limiting the number of cores to use
    if (
        worker_type_param != "manual"
        and worker_id_param
        and worker_id_param.find("_") != -1
    ):
        # ex)
        # total_cpu_num = 32, num_workers_per_node_param = 4
        # split_cpu_num = 8
        # worker_number - 1 == 0 --> [0, 1, 2, 3, 4, 5, 6, 7]
        # worker_number - 1 == 1 --> [8, 9, 10, 11, 12, 13, 14, 15]
        proc = psutil.Process(PARENT_PROCESS_ID)
        try:
            # Use the appended worker id number as worker_number
            # ex) 5wZwyCM8rxgNtERsU8znJU_1 --> extract "1" --> worker number
            worker_number = int(worker_id_param.split("_")[-1]) - 1
        except ValueError:
            logger.exception(
                "Not an expected worker ID. Cancelling CPU affinity setting"
            )
        else:
            # Note: may need to use num_cpus_to_request outside LBL
            total_cpu_num = psutil.cpu_count()
            logger.info(f"Total number of cores available: {total_cpu_num}")
            split_cpu_num = int(total_cpu_num / num_workers_per_node)
            cpu_affinity_list = list(
                range(
                    worker_number * split_cpu_num, ((worker_number + 1) * split_cpu_num)
                )
            )
            logger.info(f"Set CPU affinity to use: {cpu_affinity_list}")
            try:
                proc.cpu_affinity(cpu_affinity_list)
            except Exception as e:
                logger.exception(f"Failed to set the CPU usage limit: {e}")
                sys.exit(1)

    # Set memory upper limit
    # Todo: May need to use all free_memory on Cori and Lbl
    system_free_mem_bytes = get_free_memory()
    logger.info(
        "Total available memory (MBytes): %d"
        % (system_free_mem_bytes / 1024.0 / 1024.0)
    )

    # This available memory validation needs to executed on a compute node
    # not on a MOM node.
    if (
        worker_type_param != "manual"
        and num_workers_per_node > 1
        and not charging_account_param
    ):
        try:
            mem_per_node_to_request_byte = (
                int(mem_per_node_to_request.lower().replace("gb", "").replace("g", ""))
                * 1024.0
                * 1024.0
                * 1024.0
            )
            logger.info(
                "Requested memory for this worker (MBytes): %d"
                % (mem_per_node_to_request_byte / 1024.0 / 1024.0)
            )

            # if requested mempernode is larger than system avaiable mem space
            if system_free_mem_bytes < mem_per_node_to_request_byte:
                logger.critical("Requested memory space is not available")
                logger.critical(
                    "Available space: %d (MBytes)"
                    % (system_free_mem_bytes / 1024.0 / 1024.0)
                )
                logger.critical(
                    "Requested space: %d (MBytes)"
                    % (mem_per_node_to_request_byte / 1024.0 / 1024.0)
                )
                # Option 1
                # mem_per_node_to_request_byte = system_free_mem_bytes
                # Option 2
                raise MemoryError

            mem_limit_per_worker_bytes = int(
                mem_per_node_to_request_byte / num_workers_per_node
            )
        except Exception as e:
            logger.exception(
                f"Failed to compute the memory limit: {mem_per_node_to_request}"
            )
            logger.exception(e)
            sys.exit(1)

        try:
            soft, hard = resource.getrlimit(resource.RLIMIT_AS)
            resource.setrlimit(resource.RLIMIT_AS, (mem_limit_per_worker_bytes, hard))
            logger.info(
                "Set the memory usage upper limit (MBytes): %d"
                % (mem_limit_per_worker_bytes / 1024.0 / 1024.0)
            )
        except Exception as e:
            logger.exception(
                f"Failed to set the memory usage limit: {mem_per_node_to_request}"
            )
            logger.exception(e)
            sys.exit(1)

    job_time_to_request = (
        job_time_to_request_param
        if job_time_to_request_param
        else CONFIG.configparser.get("SLURM", "jobtime")
    )
    constraint = (
        constraint_param
        if constraint_param
        else CONFIG.configparser.get("SLURM", "constraint")
    )
    charging_account = (
        charging_account_param
        if charging_account_param
        else CONFIG.configparser.get("SLURM", "charge_accnt")
    )
    qos = qos_param if qos_param else CONFIG.configparser.get("SLURM", "qos")
    partition = (
        partition_param
        if partition_param
        else CONFIG.configparser.get("SLURM", "partition")
    )

    global THIS_WORKER_TYPE
    THIS_WORKER_TYPE = worker_type_param
    job_name = "jtm_worker_" + pool_name_param

    # Set task queue name
    if heartbeat_interval_param:
        hearbeat_interval = heartbeat_interval_param

    # Start hb receive thread
    tp_name = (
        pool_name_param
        if pool_name_param
        else CONFIG.configparser.get("JTM", "pool_name")
    )

    assert pool_name_param is not None, "User pool name is not set"
    inner_task_request_queue = (
        CONFIG.configparser.get("JTM", "jtm_inner_request_q") + "." + pool_name_param
    )

    worker_clone_time_rate = (
        worker_clone_time_rate_param
        if worker_clone_time_rate_param
        else CONFIG.configparser.getfloat("JTM", "clone_time_rate")
    )

    if THIS_WORKER_TYPE in ("dynamic"):
        assert (
            cluster_name_param != "" and cluster_name_param != "local"
        ), "Static or dynamic worker needs a cluster setting (-cl)."

    slurm_job_id = slurm_job_id_param
    cluster_name = None
    if cluster_name_param:
        cluster_name = cluster_name_param.lower()

    if (
        cluster_name == "cori"
        and mem_per_cpu_to_request != ""
        and float(
            mem_per_cpu_to_request.replace("GB", "").replace("G", "").replace("gb", "")
        )
        > 1.0
    ):
        logger.critical(
            "--mem-per-cpu in Cori shouldn't be larger than 1GB. User '--mem' instead."
        )
        sys.exit(1)

    logger.info("Task queue name: %s", inner_task_request_queue)
    logger.info("Worker type: %s", THIS_WORKER_TYPE)

    env_act = CONFIG.configparser.get("JTM", "env_activation")

    if slurm_job_id == 0 and THIS_WORKER_TYPE in ["dynamic"]:
        batch_job_script_file = os.path.join(
            job_script_dir_name,
            "jtm_%s_worker_%s.job" % (THIS_WORKER_TYPE, UNIQ_WORKER_ID),
        )
        batch_job_script_str = ""
        batch_job_misc_params = ""

        worker_config = CONFIG.config_file if CONFIG else ""
        if CONFIG.configparser.get("JTM", "worker_config_file"):
            worker_config = CONFIG.configparser.get("JTM", "worker_config_file")

        if cluster_name in ("cori", "jgi", "tahoma"):

            with open(batch_job_script_file, "w") as jf:
                batch_job_script_str += "#!/bin/bash -l"

                if cluster_name in ("cori"):
                    if num_nodes_to_request_param:
                        batch_job_script_str += f"""
#SBATCH -N {num_nodes_to_request}
#SBATCH --mem={mem_per_node_to_request}"""

                        batch_job_misc_params += (
                            f" -N {num_nodes_to_request} -m {mem_per_node_to_request} "
                        )

                        if num_cores_to_request_param:
                            batch_job_script_str += f"""
#SBATCH -c {num_cpus_to_request}"""
                            batch_job_misc_params += f" -c {num_cpus_to_request} "

                    else:
                        batch_job_script_str += f"""
#SBATCH -c {num_cpus_to_request}"""
                        batch_job_misc_params += f" -c {num_cpus_to_request} "

                        if mem_per_node_to_request:
                            batch_job_script_str += f"""
#SBATCH --mem={mem_per_node_to_request}"""
                            batch_job_misc_params += f" -m {mem_per_node_to_request} "
                        else:
                            batch_job_script_str += f"""
#SBATCH --mem-per-cpu={mem_per_cpu_to_request}"""
                            batch_job_misc_params += f" -mc {mem_per_cpu_to_request} "

                        if worker_id_param:
                            batch_job_misc_params += " -wi %(worker_id)s_${i} " % dict(
                                worker_id=UNIQ_WORKER_ID
                            )

                    ###########################
                    if 1:
                        # Need to set both --qos=genepool (or genepool_shared) _and_ -A fungalp
                        # OR
                        # no qos _and_ -A m342 _and_ -C haswell

                        # Note: currently constraint in ["haswell" | "knl"]
                        if constraint == "haswell":
                            if qos_param:
                                batch_job_script_str += f"""
#SBATCH -q {qos}"""
                                batch_job_misc_params += f" -q {qos} "
                            else:
                                batch_job_script_str += f"""
#SBATCH -q {qos}"""
                            batch_job_script_str += """
#SBATCH -C haswell"""

                            if charging_account == "m342":
                                batch_job_misc_params += " -A m342 "

                            batch_job_script_str += f"""
#SBATCH -A {charging_account}"""

                        elif constraint == "knl":
                            # Note: Basic KNL setting = "-q regular -A m342 -C knl"
                            #
                            # Note: KNL MCDRAM setting -> cache or flat
                            #  cache mode - MCDRAM is configured entirely as a last-level cache (L3)
                            #  flat mode - MCDRAM is configured entirely as addressable memory
                            #  ex) #SBATCH -C knl,quad,cache
                            #  ex) #SBATCH -C knl,quad,flat
                            #      --> srun <srun options> numactl -p 1 yourapplication.x
                            #
                            # Note: for knl, we should use m342
                            #
                            # Note: for knl, charging_account can be set via runtime (like lanl, m3408)
                            #
                            batch_job_script_str += f"""
#SBATCH -C knl
#SBATCH -A {charging_account}
#SBATCH -q {qos}"""
                            batch_job_misc_params += f" -A {charging_account} -q {qos} "

                        elif constraint == "skylake":
                            # Example usage with skylakte for Brian F.
                            # 120G
                            # ======================
                            # -t 48:00:00 -c 16 --job-name=mga-627530 --mem=115G --qos=genepool_special
                            # --exclusive -A gtrqc
                            #
                            # 250G
                            # ======================
                            # -t 96:00:00 -c 72 --job-name=mga-627834 --mem=240G -C skylake --qos=jgi_exvivo
                            # -A gtrqc
                            #
                            # 500G
                            # ======================
                            # -t 96:00:00 -c 72 --job-name=mga-627834 --mem=240G -C skylake --qos=jgi_exvivo
                            # -A gtrqc
                            #
                            # SBATCH examples
                            # 500G
                            # module load esslurm; sbatch --qos=jgi_exvivo -C skylake -A pkasmb -N 1
                            # -t 48:00:00 -D $PWD --wrap=''
                            #
                            # 250G
                            # module load esslurm; sbatch --qos=jgi_shared --mem=250G --cpus-per-task=12
                            # -C skylake -A pkasmb
                            # -N 1 -t 12:00:00 -D $PWD --wrap=""

                            batch_job_script_str += f"""
#SBATCH -N {num_nodes_to_request}
#SBATCH -C skylake
#SBATCH -A {charging_account}
#SBATCH -q {qos}"""
                            batch_job_misc_params += f" -A {charging_account} -q {qos}"

                        excl_param = ""
                        if constraint != "skylake":
                            excl_param = "#SBATCH --exclusive"

                        tq_param = ""
                        if pool_name_param:
                            tq_param = f"-p {pool_name_param}"

                        batch_job_script_str += """
#SBATCH -t %(wall_time)s
#SBATCH --job-name=%(job_name)s
#SBATCH -o %(job_dir)s/jtm_%(worker_type)s_worker_%(worker_id)s.out
#SBATCH -e %(job_dir)s/jtm_%(worker_type)s_worker_%(worker_id)s.err
%(exclusive)s

module unload python
%(env_activation_cmd)s
%(export_jtm_config_file)s

%(pagurus)s \
--move \
--user $USER \
--path %(pmetrics_path)s \
--outfile %(job_name)s_$SLURM_JOB_ID.csv &

PID=$!
sleep 2

for i in {1..%(num_workers_per_node)d}
do
    echo "jobid: $SLURM_JOB_ID"
    jtm %(set_jtm_config_file)s %(debug)s worker --slurm_job_id $SLURM_JOB_ID \
-cl %(nersc_cluster_name)s \
-wt %(worker_type)s \
-t %(wall_time)s \
--clone_time_rate %(clone_time_rate)f %(task_queue)s \
--num_worker_per_node %(num_workers_per_node)d \
-C %(constraint)s \
-m %(mem)s \
%(other_params)s &
sleep 0.5
done
wait
""" % dict(
                            debug="--debug" if DEBUG else "",
                            wall_time=job_time_to_request,
                            job_dir=job_script_dir_name,
                            worker_id=UNIQ_WORKER_ID,
                            worker_type=THIS_WORKER_TYPE,
                            clone_time_rate=worker_clone_time_rate,
                            task_queue=tq_param,
                            num_workers_per_node=num_workers_per_node,
                            env_activation_cmd=env_act,
                            other_params=batch_job_misc_params,
                            constraint=constraint,
                            mem=mem_per_node_to_request,
                            nersc_cluster_name=cluster_name,
                            job_name=job_name,
                            exclusive=excl_param,
                            export_jtm_config_file="export JTM_CONFIG_FILE=%s" % worker_config,
                            set_jtm_config_file="--config=%s" % worker_config,
                            pagurus=CONFIG.configparser.get("PERFORMANCE_METRICS", "script"),
                            pmetrics_path=CONFIG.configparser.get("PERFORMANCE_METRICS", "output_dir"),
                            pmetrics_log=CONFIG.configparser.get("PERFORMANCE_METRICS", "log_file"),
                        )

                elif cluster_name in ("jgi"):
                    if worker_id_param:
                        batch_job_misc_params += " -wi %(worker_id)s_${i} " % dict(
                            worker_id=UNIQ_WORKER_ID
                        )

                    if num_cores_to_request_param:
                        batch_job_script_str += f"""
#SBATCH --cpus-per-task {num_cpus_to_request}"""
                        batch_job_misc_params += f" -c {num_cpus_to_request}"

                    tp_param = ""
                    if pool_name_param:
                        tp_param = f"-p {pool_name_param}"

                    mnode_param = f"#SBATCH --mem={mem_per_node_to_request}"

                    batch_job_script_str += """
#SBATCH --time=%(wall_time)s
#SBATCH --job-name=%(job_name)s
#SBATCH --partition=%(partition_name)s
#SBATCH --qos=%(qosname)s
#SBATCH --account=%(charging_account)s
#SBATCH --nodes=%(num_nodes_to_request)d
%(mem_per_node_setting)s
#SBATCH -o %(job_dir)s/jtm_%(worker_type)s_worker_%(worker_id)s.out
#SBATCH -e %(job_dir)s/jtm_%(worker_type)s_worker_%(worker_id)s.err

%(env_activation_cmd)s
%(export_jtm_config_file)s
for i in {1..%(num_workers_per_node)d}
do
    echo "jobid: $SLURM_JOB_ID"
    jtm %(set_jtm_config_file)s %(debug)s worker --slurm_job_id $SLURM_JOB_ID \
-cl %(lblit_cluster_name)s \
-wt %(worker_type)s \
-t %(wall_time)s \
--clone_time_rate %(clone_time_rate)f %(task_queue)s \
--num_worker_per_node %(num_workers_per_node)d \
-m %(mem)s \
%(other_params)s &
sleep 0.5
done
wait
""" % dict(
                        debug="--debug" if DEBUG else "",
                        wall_time=job_time_to_request,
                        job_name=job_name,
                        partition_name=partition,
                        qosname=qos,
                        charging_account=charging_account,
                        num_nodes_to_request=num_nodes_to_request,
                        mem_per_node_setting=mnode_param,
                        worker_id=UNIQ_WORKER_ID,
                        job_dir=job_script_dir_name,
                        env_activation_cmd=env_act,
                        num_workers_per_node=num_workers_per_node,
                        mem=mem_per_node_to_request,
                        lblit_cluster_name=cluster_name,
                        worker_type=THIS_WORKER_TYPE,
                        clone_time_rate=worker_clone_time_rate,
                        task_queue=tp_param,
                        other_params=batch_job_misc_params,
                        export_jtm_config_file="export JTM_CONFIG_FILE=%s"
                        % worker_config,
                        set_jtm_config_file="--config=%s" % worker_config,
                    )

                elif cluster_name in ("tahoma"):
                    if worker_id_param:
                        batch_job_misc_params += " -wi %(worker_id)s_${i} " % dict(
                            worker_id=UNIQ_WORKER_ID
                        )

                    tp_param = ""
                    if pool_name_param:
                        tp_param = f"-p {pool_name_param}"

                    if partition_param:
                        batch_job_script_str += f"""
#SBATCH --partition={partition_param}"""

                    mnode_param = f"#SBATCH --mem={mem_per_node_to_request}"

                    batch_job_script_str += """
#SBATCH --account=%(charging_account)s
#SBATCH --nodes=%(num_nodes_to_request)d
#SBATCH --ntasks-per-node 16
#SBATCH --time=%(wall_time)s
#SBATCH --job-name=%(job_name)s
%(mem_per_node_setting)s
#SBATCH -o %(job_dir)s/jtm_%(worker_type)s_worker_%(worker_id)s.out
#SBATCH -e %(job_dir)s/jtm_%(worker_type)s_worker_%(worker_id)s.err

%(env_activation_cmd)s
%(export_jtm_config_file)s
for i in {1..%(num_workers_per_node)d}
do
    echo "jobid: $SLURM_JOB_ID"
    jtm %(set_jtm_config_file)s %(debug)s worker --slurm_job_id $SLURM_JOB_ID \
-cl %(emsl_cluster_name)s \
-wt %(worker_type)s \
-t %(wall_time)s \
--clone_time_rate %(clone_time_rate)f %(task_queue)s \
--num_worker_per_node %(num_workers_per_node)d \
-m %(mem)s \
%(other_params)s &
sleep 0.5
done
wait
""" % dict(
                        debug="--debug" if DEBUG else "",
                        wall_time=job_time_to_request,
                        job_name=job_name,
                        charging_account=charging_account,
                        num_nodes_to_request=num_nodes_to_request,
                        worker_id=UNIQ_WORKER_ID,
                        job_dir=job_script_dir_name,
                        env_activation_cmd=env_act,
                        num_workers_per_node=num_workers_per_node,
                        mem_per_node_setting=mnode_param,
                        mem=mem_per_node_to_request,
                        emsl_cluster_name=cluster_name,
                        worker_type=THIS_WORKER_TYPE,
                        clone_time_rate=worker_clone_time_rate,
                        task_queue=tp_param,
                        other_params=batch_job_misc_params,
                        export_jtm_config_file="export JTM_CONFIG_FILE=%s"
                        % worker_config,
                        set_jtm_config_file="--config=%s" % worker_config,
                    )

                jf.writelines(batch_job_script_str)

            os.chmod(batch_job_script_file, 0o775)

            if dry_run:
                print(batch_job_script_str)
                sys.exit(0)

            sbatch_cmd = f"sbatch --parsable {batch_job_script_file}"

            if constraint == "skylake" or qos in ("jgi_exvivo", "jgi_shared"):
                sbatch_cmd = "module load esslurm; " + sbatch_cmd
                logger.debug(f"skylake sbatch: {sbatch_cmd}")
            _, _, ec = run_sh_command(sbatch_cmd, log=logger)
            return ec

        elif cluster_name == "aws":
            pass

    logger.info("Waiting for a request...")
    logger.debug(f"Main pid = {PARENT_PROCESS_ID}")

    pid_list = []

    def signal_handler(signum, frame):
        proc_clean_exit(pid_list)

    signal.signal(signal.SIGTERM, signal_handler)

    # Start task termination proc
    max_retries = CONFIG.configparser.getint("JTM", "max_retries")
    try:
        task_kill_proc_hdl = mp.Process(
            target=TaskTerminator(config=CONFIG, max_retries=max_retries).start
        )
        task_kill_proc_hdl.start()
        pid_list.append(task_kill_proc_hdl)
    except Exception as e:
        logger.exception(f"recv_task_kill_request_proc: {e}")
        proc_clean_exit(pid_list)
        raise

    # Start send_hb_to_client_proc proc
    try:
        send_hb_to_client_proc_hdl = mp.Process(
            target=send_hb_to_client_proc,
            args=(
                hearbeat_interval,
                slurm_job_id,
                mem_per_node_to_request,
                mem_per_cpu_to_request,
                num_cpus_to_request,
                job_time_to_request,
                tp_name,
                num_workers_per_node,
                show_resource_log,
            ),
        )
        send_hb_to_client_proc_hdl.start()
        pid_list.append(send_hb_to_client_proc_hdl)
    except Exception as e:
        logger.exception(f"send_hb_to_client_proc: {e}")
        proc_clean_exit(pid_list)
        raise

    logger.info(
        "Start sending my heartbeat to the client in every %d sec to %s"
        % (hearbeat_interval, CONFIG.configparser.get("JTM", "worker_hb_q_postfix"))
    )

    # Start task runner proc
    try:
        process_task_proc_hdl = mp.Process(
            target=TaskRunner(config=CONFIG, max_retries=max_retries).start,
            args=(inner_task_request_queue,),
        )
        process_task_proc_hdl.start()
        pid_list.append(process_task_proc_hdl)
    except Exception as e:
        logger.exception(f"process_task_proc: {e}")
        proc_clean_exit(pid_list)
        raise

    return 0
