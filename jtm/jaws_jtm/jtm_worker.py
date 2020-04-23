#! /usr/bin/env python
# pylint: disable=C0111,C0103,R0205
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)

"""

jtm worker


Example scenario
1. jtm_submit sends a msg to "jtm_main_exchange" with "jtm_task_request_queue" tag.
2. jtm listens to "jtm_task_request_queue" which is bound to "jtm_main_exchange"
3. when a task is ready, jtm takes it and sends it to one of the workers
   (to jgi_jtm_inner_main_exchange)
4. workers listen to jgi_jtm_inner_request_queue which is bound to "jgi_jtm_inner_main_exchange"
5. when a task is done, a worker sends a result msg to "jgi_jtm_inner_main_exchange" with
   "jgi_jtm_inner_result_queue" tag
6. jtm listens to "jgi_jtm_inner_result_queue" queue; when a result is ready, takes and updates
   tables

"""
import os
import sys
import multiprocessing as mp
import shortuuid
import datetime
import json
import threading
import subprocess
import socket
import time
import pika
import resource
import psutil
import functools
import signal

from jaws_jtm.common import setup_custom_logger, logger
from jaws_jtm.lib.rabbitmqconnection import RmqConnectionHB
from jaws_jtm.lib.resourceusage import get_cpu_load, \
    get_runtime, get_pid_tree, get_virtual_memory_usage, \
    get_resident_memory_usage, get_total_mem_usage_per_node, \
    get_num_workers_on_node, get_free_memory
from jaws_jtm.lib.run import make_dir, run_sh_command
from jaws_jtm.lib.msgcompress import zdumps, zloads


# This ipc pipe is to send task_id (by do_work())to send_hb_to_client_proc()
# when a task is requested
PIPE_TASK_ID_SEND, PIPE_TASK_ID_RECV = mp.Pipe()
DEBUG = False

# -------------------------------------------------------------------------------
# Globals
# -------------------------------------------------------------------------------
# To share user task proc id
USER_PROC_PROC_ID = mp.Value("i", 0)
# To share this worker life left
WORKER_LIFE_LEFT_IN_MINUTE = mp.Value("i", 0)
ALLOW_EMPTY_OUTPUT_FILE = False  # False ==> if size(output file)==0, return error code
UNIQ_WORKER_ID = str(shortuuid.uuid())
IS_CLIENT_ALIVE = False
WORKER_START_TIME = datetime.datetime.now()
PARENT_PROCESS_ID = os.getpid()  # parent process id
THIS_WORKER_TYPE = None


# -------------------------------------------------------------------------------
def run_user_task(msg_unzipped, return_msg, ch):
    """
    Run a user command in msg_zipped_to_send
    :param msg_unzipped: uncompressed msg from client
    :param return_msg: msg to return
    :param ch: rmq channel
    :return:
    """
    # Uncompress msg to get a taska
    task_id = msg_unzipped["task_id"]
    user_task_cmd = msg_unzipped["user_cmd"]
    out_files = msg_unzipped["output_files"]
    task_type = msg_unzipped["task_type"]

    return_msg["task_id"] = task_id
    return_msg["user_cmd"] = user_task_cmd
    return_msg["task_type"] = task_type

    if "stdout" in msg_unzipped:
        user_task_cmd + " > %s" % (msg_unzipped["stdout"])
    if "stderr" in msg_unzipped:
        user_task_cmd + " 2>%s" % (msg_unzipped["stderr"])

    # Run the task
    logger.info("Running task ID %d...", task_id)

    p = None
    time_out_in_minute = 0

    if THIS_WORKER_TYPE != "manual":
        # wait until WORKER_LIFE_LEFT_IN_MINUTE is updated
        if WORKER_LIFE_LEFT_IN_MINUTE.value <= 0:
            while True:
                logger.debug("worker life: %d", WORKER_LIFE_LEFT_IN_MINUTE.value)
                if WORKER_LIFE_LEFT_IN_MINUTE.value > 0:
                    break
                ch._connection.process_data_events(time_limit=1.0)

        # ex) WORKER_LIFE_LEFT_IN_MINUTE = 20min and TASK_KILL_TIMEOUT_MINUTE = 10min
        # timeout will be set as 10min
        # TASK_KILL_TIMEOUT_MINUTE is a extra housekeeping time after explicitly
        # terminate a task.
        logger.debug("worker life: %d", WORKER_LIFE_LEFT_IN_MINUTE.value)
        time_out_in_minute = int(WORKER_LIFE_LEFT_IN_MINUTE.value - CONFIG.constants.TASK_KILL_TIMEOUT_MINUTE)
        logger.info("Timeout in minute: %d", time_out_in_minute)

    proc_return_code = -1
    done_flags = CONFIG.constants.DONE_FLAGS
    logger.debug("Start subprocess to run a task.")
    try:
        p = subprocess.Popen(user_task_cmd.split(),
                             env=os.environ,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
    except MemoryError:
        logger.exception("Exception: Out of memory, %s", user_task_cmd)
        proc_return_code = done_flags["failed with out-of-mem"]
    except FileNotFoundError:
        logger.exception("Exception: Input file or command not found, %s", user_task_cmd)
        proc_return_code = 2
    except Exception as detail:
        logger.exception("Exception: Failed to run user command, %s", user_task_cmd)
        logger.exception("Detail: %s", str(detail))
    else:
        # Set USER_PROC_PROC_ID = forked child process id (this value will be sent
        # to send_hb_to_client function)
        USER_PROC_PROC_ID.value = p.pid
        logger.debug("USER_PROC_PROC_ID.value %d" % USER_PROC_PROC_ID.value)
        logger.debug("User command: %s", user_task_cmd)

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
            #
            proc_return_code = done_flags["failed with timeout"]
        else:
            p.wait()
            p.poll()
            proc_return_code = p.returncode
    else:
        logger.error("subprocess call failed")

    # Prepare result to send back
    return_msg["out_files"] = out_files
    return_msg["worker_id"] = UNIQ_WORKER_ID
    return_msg["host_name"] = socket.gethostname()

    if USER_PROC_PROC_ID.value == -9:
        return_msg["done_flag"] = done_flags["failed with user termination"]
        return_msg["ret_msg"] = "Task cancelled."
    else:
        if proc_return_code == done_flags["failed with timeout"]:
            logger.info("User task timeout. Not enough time left for the worker.")
            return_msg["done_flag"] = str(done_flags["failed with timeout"])
            return_msg["ret_msg"] = "User task timeout"
        elif proc_return_code == done_flags["failed with out-of-mem"]:
            logger.info("User task out-of-mem.")
            return_msg["done_flag"] = str(done_flags["failed with out-of-mem"])
            return_msg["ret_msg"] = "User task out-of-mem"
        elif proc_return_code == 2:  # system code
            logger.info("input file or command not found.")
            return_msg["done_flag"] = done_flags["failed with input file or command not found"]
            return_msg["ret_msg"] = "Input file or command not found."
        elif proc_return_code == 0:  # system code
            logger.info("Task# %s completed!" % (str(task_id)))

            # Output file checking
            if out_files:
                ofs = out_files.split(",")
                logger.debug("Number of output files = %d.", len(ofs))
                out_file_list = []

                for i in range(len(ofs)):
                    out_file_list.append(ofs[i])

                ret, file_size = check_output(out_file_list,
                                              CONFIG.configparser.getfloat("JTM", "file_check_interval"),
                                              CONFIG.configparser.getint("JTM", "file_checking_max_trial"),
                                              CONFIG.configparser.getfloat("JTM", "file_check_int_inc"))

                if not ret:
                    ret_msg_str = "Failed to check output file(s): %s, file size = %s." % (ofs, file_size)
                    logger.critical(ret_msg_str)
                    return_msg["done_flag"] = done_flags["failed to check output file(s)"]
                    return_msg["ret_msg"] = ret_msg_str
                else:
                    return_msg["done_flag"] = done_flags["success with correct output file(s)"]
                    return_msg["ret_msg"] = "Output file checking is OK."
            else:
                return_msg["done_flag"] = done_flags["success"]
                return_msg["ret_msg"] = "No file(s) to check."
        else:
            logger.critical("Failed to execute a task, %s. Non-zero exit code. stdout = %s."
                            % (user_task_cmd, stdout_str))
            return_msg["done_flag"] = done_flags["failed to run user command"]
            return_msg["ret_msg"] = stdout_str

    logger.info("Reply msg prepared with result: %s" % str(return_msg))


# -------------------------------------------------------------------------------
def check_output(out_files, out_file_check_wait_time=3,
                 max_trial=3, out_file_check_wait_time_increase=1.5):
    """
    Check 1) existence, 2) size>0 for each file in out_files
    :param out_files: list of absolute paths to the output files to check
    :param out_file_check_wait_time: sleep time between output file checking before retiral
    :param max_trial: max trial for checking
    :param out_file_check_wait_time_increase: wait time increase for retrial
    :return:
    """
    file_size = 0
    b_is_file_found = False
    trial = 1

    for a_file in out_files:
        logger.info("Output file check: %s", a_file)
        logger.debug("Output file check: %s", os.path.expandvars(a_file))
        a_file = os.path.expandvars(a_file)

        file_size = 0
        b_is_file_found = False

        while trial < max_trial:
            logger.info("Output file checking. Trial# = %d", trial)

            # First, check file existence
            # os.path.exists returns if it is a valid path(check for directory or file, both)
            # and os.path.isfile(checks for only file, not directory) returns if it is a file
            b_is_file_exist = os.path.isfile(a_file)

            # If exist, check file size
            if b_is_file_exist:
                file_size = os.path.getsize(a_file)
                if file_size == 0:
                    logger.warning("File, %s is zero size.", a_file)
                if file_size > 0:
                    b_is_file_found = True
                    logger.info("Output file '%s' is OK.", a_file)
                    break
            else:
                logger.info("Outout file not found.")

            # Wait for initial wait time
            time.sleep(out_file_check_wait_time)

            # Increase the wait time
            out_file_check_wait_time *= out_file_check_wait_time_increase
            trial += 1

    return b_is_file_found, file_size


# -------------------------------------------------------------------------------
def send_hb_to_client_proc(interval, slurm_job_id, mem_per_node, mem_per_core,
                           num_cores, job_time, task_queue_name, pool_name,
                           nwpn, exch_name, worker_hb_queue):
    """
    Send heartbeats to the client
    :param interval: time interval to send heartbeats to the client
    :param slurm_job_id: SLURM job id
    :param mem_per_node: memory request per node
    :param mem_per_core: memory request per core
    :param num_cores: number of cores
    :param job_time: wallclocktime
    :param task_queue_name: task queue name
    :param pool_name: pool name
    :param nwpn: number of workers per node
    :param exch_name
    :param worker_hb_queue: worker hb queue postfix
    :return:
    """
    # Remote broker (mq.nersc.gov)
    # rmq_conn = RmqConnection()
    # with heartbeat_interval=0
    # ref) http://stackoverflow.com/questions/14572020/handling-long-running-tasks-in-pika-rabbitmq,
    # http://stackoverflow.com/questions/34721178/pika-blockingconnection-rabbitmq-connection-closed
    rmq_conn = RmqConnectionHB(config=CONFIG)
    conn = rmq_conn.open()
    ch = conn.channel()
    host_name = socket.gethostname()

    ch.exchange_declare(exchange=exch_name,
                        exchange_type="direct",
                        passive=False,
                        durable=False,
                        auto_delete=False)

    task_id = 0
    proc_id_list_merged = []

    while True:
        try:
            # Todo: make it optional to send the resource info on hb message
            #
            # To check out-of-mem
            # 1. On genepool node
            #    Need to get ram.c and v_mem.c, and compare them with sum(mem usages
            #    from all processes) to check out-of-mem
            # 2. On Mendel node
            #    Use "free" call for getting mem usage (%) for the node
            #    If >90% used, kill a process (selection strategy is needed)
            #
            # NOTE: must cope with the fast mem consumption at the begining of the process

            if USER_PROC_PROC_ID.value == 0:
                child_pid = PARENT_PROCESS_ID
            else:
                child_pid = int(USER_PROC_PROC_ID.value)

            # Collect pids from process tree
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

            vmem_usage_list = []
            rmem_usage_list = []

            try:
                vmem_usage_list.extend([get_virtual_memory_usage(pid, 0.0, False) for pid in proc_id_list_merged])
            except ValueError:
                logger.exception("ValueError: Failed to collect VM memory usage.")
            except UnboundLocalError:
                logger.exception("UnboundLocalError: No entry in process id list.")

            try:
                rmem_usage_list.extend([get_resident_memory_usage(pid, 0.0, False) for pid in proc_id_list_merged])
            except ValueError:
                logger.exception("ValueError: Failed to collect RES memory usage.")
            except UnboundLocalError:
                logger.exception("UnboundLocalError: No entry in process id list.")

            # Collect cpu_usages for all pids in the tree and get max()
            try:
                cpu_load_list = [get_cpu_load(pid) for pid in proc_id_list_merged]
            except Exception as e:
                logger.exception("get_cpu_load() exception: {}".format(e))
                raise

            max_cpu_load = max(cpu_load_list) if len(cpu_load_list) > 0 else 0.0

            # Collect mem_usages for all pids in the tree and get sum()
            rmem_usage = "%.1f" % sum(rmem_usage_list)
            vmem_usage = "%.1f" % sum(vmem_usage_list)

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
            perc_used_mem = "%.1f" % get_total_mem_usage_per_node()
            num_workers_on_node = get_num_workers_on_node(CONFIG)

            # Check if there is any task id in the ipc pipe
            if PIPE_TASK_ID_RECV.poll():
                task_id = PIPE_TASK_ID_RECV.recv()

            today = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            end_date_time = None

            try:
                # Note: this only works with SLURM. For other HPCs/Clouds, methods to get the
                #  remaining life are needed! ==> resolved!
                if slurm_job_id:
                    # NEW
                    # jobtime2 = datetime.datetime.strptime(job_time, '%H:%M:%S')
                    # end_date_time = WORKER_START_TIME + timedelta(seconds=jobtime2.second+jobtime2.minute*60+
                    # jobtime2.hour*3600)
                    # 24:00:00 --> seconds
                    job_runtime_in_sec = int(job_time.split(':')[0]) * 3600 + \
                                         int(job_time.split(':')[1]) * 60 + \
                                         int(job_time.split(':')[2])
                    end_date_time = WORKER_START_TIME + datetime.timedelta(seconds=int(job_runtime_in_sec))
                    delta = end_date_time - datetime.datetime.now()
                    WORKER_LIFE_LEFT_IN_MINUTE.value = int(divmod(delta.total_seconds(), 60)[0])
            except Exception as e:
                logger.critical("Something wrong in computing remaining wall clock time: %s", e)
                raise

            # Get the total number of tasks in the task queue for this pool
            hb_queue_len = ch.queue_declare(queue=task_queue_name,
                                            durable=True,
                                            exclusive=False,
                                            auto_delete=True).method.message_count

            # To decrease the hb message size, changed to int key value
            try:
                ip_address = socket.gethostbyname(host_name)
            except Exception:
                ip_address = None

            ncore_param = mp.cpu_count()
            worker_t = CONFIG.constants.WORKER_TYPE
            if worker_t[THIS_WORKER_TYPE] > 0:
                ncore_param = num_cores

            hb_msg = CONFIG.constants.HB_MSG
            msg_dict_to_send = {hb_msg["child_pid"]: child_pid,
                                hb_msg["clone_time_rate"]: 0.0,
                                hb_msg["cpu_load"]: max_cpu_load,
                                hb_msg["end_date"]: today,  # Note: for dynamic worker endDate update
                                hb_msg["host_name"]: host_name,
                                hb_msg["ip_address"]: ip_address,
                                hb_msg["job_time"]: job_time if worker_t[THIS_WORKER_TYPE] > 0 else None,
                                hb_msg["jtm_host_name"]: CONFIG.configparser.get("SITE", "jtm_host_name"),
                                hb_msg["life_left"]: WORKER_LIFE_LEFT_IN_MINUTE.value,
                                hb_msg["mem_per_core"]: mem_per_core if worker_t[THIS_WORKER_TYPE] > 0 else "",
                                hb_msg["mem_per_node"]: mem_per_node if worker_t[THIS_WORKER_TYPE] > 0 else "",
                                hb_msg["num_cores"]: ncore_param,
                                hb_msg["num_tasks"]: hb_queue_len,
                                hb_msg["num_workers_on_node"]: num_workers_on_node,
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
                                hb_msg["worker_type"]: worker_t[THIS_WORKER_TYPE],
                                hb_msg["nwpn"]: nwpn}

            msg_zipped_to_send = zdumps(json.dumps(msg_dict_to_send))
            ch.basic_publish(exchange=exch_name,
                             routing_key=worker_hb_queue,
                             body=msg_zipped_to_send)
        except Exception as e:
            logger.critical("Something wrong with send_hb_to_client(): %s", e)
            ch.close()
            conn.close()
            raise

        # Todo: Need to start with shorted interval like (0.5-1 sec) for about 30 sec
        #  from the beginning and then use the user specified interval
        conn.sleep(interval)

    # unreachable
    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def recv_task_kill_request_proc():
    """
    Wait for task termination request from JTM

    :return:
    """

    rmq_conn = RmqConnectionHB(config=CONFIG)
    conn = rmq_conn.open()
    ch = conn.channel()
    exch_name = CONFIG.configparser.get("JTM", "jtm_task_kill_exch")
    queue_name = CONFIG.configparser.get("JTM", "jtm_task_kill_q")
    routing_key = UNIQ_WORKER_ID

    ch.exchange_declare(exchange=exch_name,
                        exchange_type="fanout",
                        durable=True,
                        auto_delete=False)
    ch.queue_declare(queue=queue_name,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)
    ch.queue_bind(exchange=exch_name,
                  queue=queue_name,
                  routing_key=routing_key)

    def on_kill(ch, method, props, body):
        msg_unzipped = json.loads(zloads(body))
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
                    logger.info("Executing {} for taskID, {}".format(kill_cmd, msg_unzipped["task_id"]))
                    _, _, ec = run_sh_command(kill_cmd, log=logger)
                    if ec == 0:
                        logger.info("Successfully terminate a user task process.")
                    else:
                        logger.warning("User process not found. Ignore the termination command, %s"
                                       % (kill_cmd))
                        # Todo: Failed to terminate a user process for some reason. How to deal
                        #  with this case?
                        # ch.basic_reject(delivery_tag=method.delivery_tag, requeue=True)

                # Kill the main child process
                # Note: can consider to use "pkill -9 -P ppid" to kill the family
                kill_cmd = "kill -9 %d" % msg_unzipped["child_pid"]
                logger.info("Executing {} for taskID, {}".format(kill_cmd, msg_unzipped["task_id"]))
                _, _, ec = run_sh_command(kill_cmd, log=logger)
                if ec == 0:
                    logger.info("Successfully terminate a user task process.")
                else:
                    logger.warning("User process not found. Failed to execute the command, %s" % (kill_cmd))

            else:
                logger.warning("No valid child process id to terminate.")

            ch.basic_ack(delivery_tag=method.delivery_tag)

        else:
            # ch.basic_nack(delivery_tag=method.delivery_tag)
            # If UNIQ_WORKER_ID is not for me, reject and requeue it
            logger.info("NACK sent to the broker")
            ch.basic_reject(delivery_tag=method.delivery_tag,
                            requeue=True)

    # Waiting for a task kill
    ch.basic_qos(prefetch_count=1)
    ch.basic_consume(queue=queue_name,
                     on_message_callback=on_kill,
                     auto_ack=False)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        ch.stop_consuming()
        raise


# -------------------------------------------------------------------------------
def check_processes(pid_list):
    """
    Checking if the total number of processes from the worker is NUM_WORKER_PROCS
    if not, terminate the all the proc ids
    :param pid_list: process id list
    :return:
    """
    while True:
        logger.debug("Total Number of processes of the worker = %d" % (len(pid_list)+2))
        if len(pid_list) != CONFIG.constants.NUM_WORKER_PROCS - 2:
            raise OSError(2, 'Number of processes is wrong')
        time.sleep(CONFIG.configparser.getfloat("JTM", "num_procs_check_interval"))


# -------------------------------------------------------------------------------
def ack_message(ch, reply_to, correlation_id, delivery_tag, response):
    """

    :param ch:
    :param reply_to:
    :param correlation_id:
    :param delivery_tag:
    :param response:
    :return:
    """
    # Note that `ch` must be the same pika channel instance via which
    #  the message being ACKed was retrieved (AMQP protocol constraint).
    if ch.is_open:
        logger.debug("Send ACK to the manager")
        ch.basic_ack(delivery_tag)
        logger.debug("Send result to the manager")
        ch.queue_declare(queue=reply_to,
                         durable=True,
                         exclusive=False,
                         auto_delete=True)
        ch.basic_publish(exchange=JTM_INNER_MAIN_EXCH,
                         routing_key=reply_to,
                         properties=pika.BasicProperties(
                              delivery_mode=2,
                              correlation_id=correlation_id),
                         body=response)
    else:
        # Channel is already closed, so we can't ACK this message;
        # log and/or do something that makes sense for your app in this case.
        logger.info("Sending ACK failed. Try to send 'failed' to the manager")
        ###########################################################
        # Note: This per-message & per-worker connection might
        # affect the broker performance.
        ###########################################################
        rmq_conn3 = RmqConnectionHB(config=CONFIG)
        conn3 = rmq_conn3.open()
        ch3 = conn3.channel()
        lost_connection_response = {}
        lost_connection_response["done_flag"] = str(CONFIG.constants.DONE_FLAGS["failed with lost connection"])
        lost_connection_response["ret_msg"] = "Stream connection lost"
        json_data = json.dumps(lost_connection_response)
        response = zdumps(json_data)
        ch3.queue_declare(queue=reply_to,
                          durable=True,
                          exclusive=False,
                          auto_delete=True)
        ch3.basic_publish(exchange=JTM_INNER_MAIN_EXCH,
                          routing_key=reply_to,  # use the queue which the client created
                          properties=pika.BasicProperties(
                             delivery_mode=2,  # make message persistent
                             correlation_id=correlation_id),
                          body=response)
        ch3.basic_reject(delivery_tag=delivery_tag, requeue=True)
        ch3.close()
        conn3.close()
        raise OSError(2, 'Connection lost')


# -------------------------------------------------------------------------------
def run_task(conn, ch, delivery_tag, reply_to, correlation_id, body):
    """
    A callback function whenever a message is received.

    :param conn:
    :param ch:
    :param delivery_tag:
    :param reply_to:
    :param correlation_id:
    :param body:
    :return:
    """
    msg_unzipped = json.loads(zloads(body))
    task_id = msg_unzipped["task_id"]

    # Send taskid to send_hb_to_client_proc() so that the task_id can be shown in the hb messages
    PIPE_TASK_ID_SEND.send(task_id)

    # If client sends "terminate worker"
    if task_id == -9:
        logger.info("Received TERM signal. Terminate myself.")
        #
        # Ref) http://gavinroy.com/deeper-down-the-rabbit-hole-of-message-redeli
        # Return back the message (the TERM message) in the queue
        # so that the other workers can get it.
        logger.info("NACK sent to the broker")
        ch.basic_reject(delivery_tag=delivery_tag, requeue=True)
        raise OSError(2, 'Worker termination request')

    logger.info("Received a task, %r" % (msg_unzipped,))
    logger.debug("Return queue = %s", reply_to)

    result_dict = {}

    # OLD
    run_user_task(msg_unzipped, result_dict, ch)

    json_data = json.dumps(result_dict)
    logger.debug("Result msg sent!")
    response = zdumps(json_data)

    try:
        # Note: After sending ack, the message will be deleted from RabbitMQ
        #  If this worker crashes while running a user command, this task will
        #  be sent to other workers available

        # NEW
        cb = functools.partial(ack_message,
                               ch,
                               reply_to,
                               correlation_id,
                               delivery_tag,
                               response)
        conn.add_callback_threadsafe(cb)
    except Exception as e:
        logger.critical("Something wrong in run_task.add_callback_threadsafe(): %s", e)
        raise

    # Send taskid=0 to send_hb_to_client_proc() b/c the requested task is completed
    PIPE_TASK_ID_SEND.send(0)

    # Reset child pid
    USER_PROC_PROC_ID.value = 0


# -------------------------------------------------------------------------------
def on_task_request(ch, method_frame, _header_frame, body, args):
    """
    Threaded way to consume request from manager
    :param ch:
    :param method_frame:
    :param _header_frame:
    :param body:
    :param args:
    :return:
    """
    (conn, thrds) = args
    delivery_tag = method_frame.delivery_tag
    reply_to = _header_frame.reply_to
    correlation_id = _header_frame.correlation_id
    t = threading.Thread(target=run_task,
                         args=(conn, ch,
                               delivery_tag,
                               reply_to,
                               correlation_id,
                               body))
    t.start()
    thrds.append(t)


# -------------------------------------------------------------------------------
def proc_clean(pid_list):
    """

    :param pid_list: process handle list
    :return:
    """
    for p in pid_list:
        if p is not None and p.is_alive():
            p.terminate()
    os._exit(1)


# -------------------------------------------------------------------------------
def conn_clean(conn, ch):
    """

    :param conn:
    :param ch:
    :return:
    """
    ch.stop_consuming()
    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def worker(ctx: object, heartbeat_interval_param: int, custom_log_dir: str,
           custom_job_log_dir_name: str, pool_name_param: str,
           slurm_job_id_param: int, worker_type_param: str, cluster_name_param: str,
           num_workers_per_node_param: int,
           worker_id_param: str,
           num_cores_to_request_param: int,
           mem_per_node_to_request_param: str,
           mem_per_cpu_to_request_param: str,
           job_time_to_request_param: str) -> int:

    global CONFIG
    CONFIG = ctx.obj['config']
    global DEBUG
    DEBUG = ctx.obj['debug']
    # config file has precedence
    config_debug = CONFIG.configparser.getboolean("SITE", "debug")
    if config_debug:
        DEBUG = config_debug
    global JTM_INNER_MAIN_EXCH
    JTM_INNER_MAIN_EXCH = CONFIG.configparser.get("JTM", "jtm_inner_main_exch")
    prod_mode = False
    if CONFIG.configparser.get("JTM", "run_mode") == "prod":
        prod_mode = True

    # Job dir setting
    job_script_dir_name = os.path.join(CONFIG.configparser.get("JTM", "log_dir"), "job")
    if custom_job_log_dir_name:
        job_script_dir_name = custom_job_log_dir_name
    make_dir(job_script_dir_name)

    # Log dir setting
    log_dir_name = os.path.join(CONFIG.configparser.get("JTM", "log_dir"), "log")
    if custom_log_dir:
        log_dir_name = custom_log_dir
    make_dir(log_dir_name)

    print("JTM Worker, version: {}".format(CONFIG.constants.VERSION))

    # Set uniq worker id if worker id is provided in the params
    if worker_id_param:
        global UNIQ_WORKER_ID
        UNIQ_WORKER_ID = worker_id_param

    log_level = "info"
    if DEBUG:
        log_level = "debug"

    setup_custom_logger(log_level, log_dir_name, 1, 1, worker_id=UNIQ_WORKER_ID)

    logger.info("\n*****************\nDebug mode is %s\n*****************"
                % ("ON" if DEBUG else "OFF"))

    hearbeat_interval = CONFIG.configparser.getfloat("JTM", "worker_hb_send_interval")

    logger.info("Set jtm log file location to %s", log_dir_name)
    logger.info("Set jtm job file location to %s", job_script_dir_name)
    logger.info("RabbitMQ broker: %s", CONFIG.configparser.get("RMQ", "host"))
    logger.info("Pika version: %s", pika.__version__)
    logger.info("JTM user name: %s", CONFIG.configparser.get("SITE", "user_name"))
    logger.info("Unique worker ID: %s", UNIQ_WORKER_ID)
    logger.info("\n*****************\nRun mode is %s\n*****************"
                % ("PROD" if prod_mode else "DEV"))
    logger.info("JTM config file: %s" % (CONFIG.config_file))

    # Slurm info
    num_workers_per_node = num_workers_per_node_param \
        if num_workers_per_node_param else CONFIG.configparser.getint("JTM", "num_workers_per_node")
    assert num_workers_per_node > 0
    mem_per_cpu_to_request = mem_per_cpu_to_request_param \
        if mem_per_cpu_to_request_param else CONFIG.configparser.get("SLURM", "mempercpu")
    mem_per_node_to_request = mem_per_node_to_request_param \
        if mem_per_node_to_request_param else CONFIG.configparser.get("SLURM", "mempernode")
    assert mem_per_cpu_to_request
    assert mem_per_node_to_request
    num_cpus_to_request = num_cores_to_request_param \
        if num_cores_to_request_param else CONFIG.configparser.getint("SLURM", "ncpus")

    # Set CPU affinity for limiting the number of cores to use
    if worker_type_param != "manual" and worker_id_param and worker_id_param.find('_') != -1:
        # ex)
        # total_cpu_num = 32, num_workers_per_node_param = 4
        # split_cpu_num = 8
        # worker_number - 1 == 0 --> [0, 1, 2, 3, 4, 5, 6, 7]
        # worker_number - 1 == 1 --> [8, 9, 10, 11, 12, 13, 14, 15]
        proc = psutil.Process(PARENT_PROCESS_ID)
        try:
            # Use the appended worker id number as worker_number
            # ex) 5wZwyCM8rxgNtERsU8znJU_1 --> extract "1" --> worker number
            worker_number = int(worker_id_param.split('_')[-1]) - 1
        except ValueError:
            logger.exception("Not an expected worker ID. Cancelling CPU affinity setting")
        else:
            # Note: may need to use num_cpus_to_request outside LBL
            total_cpu_num = psutil.cpu_count()
            logger.info("Total number of cores available: {}".format(total_cpu_num))
            split_cpu_num = int(total_cpu_num / num_workers_per_node)
            cpu_affinity_list = list(range(worker_number * split_cpu_num,
                                           ((worker_number + 1) * split_cpu_num)))
            logger.info("Set CPU affinity to use: {}".format(cpu_affinity_list))
            try:
                proc.cpu_affinity(cpu_affinity_list)
            except Exception as e:
                logger.exception("Failed to set the CPU usage limit: %s" % (e))
                sys.exit(1)

    # Set memory upper limit
    # Todo: May need to use all free_memory on Cori and Lbl
    system_free_mem_bytes = get_free_memory()
    logger.info("Total available memory (MBytes): %d"
                % (system_free_mem_bytes / 1024.0 / 1024.0))

    if worker_type_param != "manual" and num_workers_per_node > 1:
        try:
            mem_per_node_to_request_byte = int(mem_per_node_to_request.lower()
                                               .replace("gb", "")
                                               .replace("g", "")) * 1024.0 * 1024.0 * 1024.0
            logger.info("Requested memory for this worker (MBytes): %d"
                        % (mem_per_node_to_request_byte / 1024.0 / 1024.0))

            # if requested mempernode is larger than system avaiable mem space
            if system_free_mem_bytes < mem_per_node_to_request_byte:
                logger.critical("Requested memory space is not available")
                logger.critical("Available space: %d (MBytes)"
                                % (system_free_mem_bytes / 1024.0 / 1024.0))
                logger.critical("Requested space: %d (MBytes)"
                                % (mem_per_node_to_request_byte / 1024.0 / 1024.0))
                # Option 1: just set the max mem available
                # mem_per_node_to_request_byte = system_free_mem_bytes
                # Option 2: report error
                raise MemoryError

            MEM_LIMIT_PER_WORKER_BYTES = int(mem_per_node_to_request_byte /
                                             num_workers_per_node)
        except Exception as e:
            logger.exception("Failed to compute the memory limit: %s", mem_per_node_to_request)
            logger.exception(e)
            sys.exit(1)

        try:
            soft, hard = resource.getrlimit(resource.RLIMIT_AS)
            resource.setrlimit(resource.RLIMIT_AS, (MEM_LIMIT_PER_WORKER_BYTES, hard))
            logger.info("Set the memory usage upper limit (MBytes): %d"
                        % (MEM_LIMIT_PER_WORKER_BYTES / 1024.0 / 1024.0))
        except Exception as e:
            logger.exception("Failed to set the memory usage limit: %s", mem_per_node_to_request)
            logger.exception(e)
            sys.exit(1)

    job_time_to_request = job_time_to_request_param \
        if job_time_to_request_param else CONFIG.configparser.get("SLURM", "jobtime")

    global THIS_WORKER_TYPE
    THIS_WORKER_TYPE = worker_type_param
    job_name = "jtm_worker_" + pool_name_param

    # Set task queue name
    if heartbeat_interval_param:
        hearbeat_interval = heartbeat_interval_param

    # Start hb receive thread
    tp_name = ""
    if pool_name_param:
        tp_name = pool_name_param

    assert pool_name_param is not None and pool_name_param != "", "User pool name is not set"
    inner_task_request_queue = CONFIG.configparser.get("JTM", "jtm_inner_request_q") + "." + pool_name_param

    if THIS_WORKER_TYPE == "dynamic":
        assert cluster_name_param != "" and \
               cluster_name_param != "local", "Static or dynamic worker needs a cluster setting (-cl)."

    slurm_job_id = slurm_job_id_param
    cluster_name = cluster_name_param

    if cluster_name == "cori" and mem_per_cpu_to_request != "" \
            and float(mem_per_cpu_to_request.replace("GB", "").replace("G", "").replace("gb", "")) > 1.0:
        logger.critical("--mem-per-cpu in Cori shouldn't be larger than 1GB. User '--mem' instead.")
        sys.exit(1)

    logger.info("Task queue name: %s", inner_task_request_queue)
    logger.info("Worker type: %s", THIS_WORKER_TYPE)

    # Remote broker (rmq.nersc.gov)
    rmq_conn = RmqConnectionHB(config=CONFIG)
    conn = rmq_conn.open()
    ch = conn.channel()
    ch.exchange_declare(exchange=JTM_INNER_MAIN_EXCH,
                        exchange_type="direct",
                        passive=False,
                        durable=True,
                        auto_delete=False)

    # Declare task receiving queue (client --> worker)
    #
    # If you have a queueu that is durable, RabbitMQ will never lose our queue.
    # If you have a queue that is exclusive, then when the channel that declared
    # the queue is closed, the queue is deleted.
    # If you have a queue that is auto-deleted, then when there are no
    # subscriptions left on that queue it will be deleted.
    #
    ch.queue_declare(queue=inner_task_request_queue,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)
    ch.queue_bind(exchange=JTM_INNER_MAIN_EXCH,
                  queue=inner_task_request_queue,
                  routing_key=inner_task_request_queue)

    logger.info("Waiting for a request...")
    pid_list = []

    # Start task termination proc
    try:
        task_kill_proc_hdl = mp.Process(target=recv_task_kill_request_proc)
        task_kill_proc_hdl.start()
        pid_list.append(task_kill_proc_hdl)
    except Exception as e:
        logger.exception("recv_task_kill_request_proc: {}".format(e))
        proc_clean(pid_list)
        conn_clean(conn, ch)
        sys.exit(1)

    # Start send_hb_to_client_proc proc
    try:
        send_hb_to_client_proc_hdl = mp.Process(target=send_hb_to_client_proc,
                                                args=(hearbeat_interval,
                                                      slurm_job_id,
                                                      mem_per_node_to_request,
                                                      mem_per_cpu_to_request,
                                                      num_cpus_to_request,
                                                      job_time_to_request,
                                                      inner_task_request_queue,
                                                      tp_name,
                                                      num_workers_per_node,
                                                      CONFIG.configparser.get("JTM", "jtm_worker_hb_exch"),
                                                      CONFIG.configparser.get("JTM", "worker_hb_q_postfix")))
        send_hb_to_client_proc_hdl.start()
        pid_list.append(send_hb_to_client_proc_hdl)
    except Exception as e:
        logger.exception("send_hb_to_client_proc: {}".format(e))
        proc_clean(pid_list)
        conn_clean(conn, ch)
        sys.exit(1)

    logger.info("Start sending my heartbeat to the client in every %d sec" % (hearbeat_interval,))

    # Checking the total number of child processes
    try:
        check_processes_hdl = mp.Process(target=check_processes,
                                         args=(pid_list,))
        check_processes_hdl.start()
        pid_list.append(check_processes_hdl)
    except Exception as e:
        logger.exception("check_processes: {}".format(e))
        proc_clean(pid_list)
        conn_clean(conn, ch)
        sys.exit(1)

    def signal_handler(signum, frame):
        proc_clean(pid_list)

    signal.signal(signal.SIGTERM, signal_handler)

    # Waiting for request
    ch.basic_qos(prefetch_count=1)

    # NEW
    # Ref) https://github.com/pika/pika/blob/1.0.1/examples/basic_consumer_threaded.py
    #      https://stackoverflow.com/questions/51752890/how-to-disable-heartbeats-with-pika-and-rabbitmq
    #      https://github.com/pika/pika/blob/master/examples/basic_consumer_threaded.py
    threads = []
    on_message_callback = functools.partial(on_task_request,
                                            args=(conn, threads))
    ch.basic_consume(queue=inner_task_request_queue,
                     on_message_callback=on_message_callback)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        proc_clean()
        conn_clean()

    # Wait for all to complete
    # Note: prefetch_count=1 ==> #thread = 1
    for thread in threads:
        thread.join()

    if conn:
        conn.close()

    return 0
