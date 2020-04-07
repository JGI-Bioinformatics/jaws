#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)

"""

jtm worker


Example scenario
1. jtm_submit sends a msg to "jgi_microservice" with "jtm_task_request_queue" tag.
2. jtm listens to "jtm_task_request_queue" which is bound to "jgi_microservice"
3. when a task is ready, jtm takes it and sends it to one of the workers
   (to jgi_jtm_inner_main_exchange)
4. workers listen to jgi_jtm_inner_request_queue which is bound to "jgi_jtm_inner_main_exchange"
5. when a task is done, a worker sends a result msg to "jgi_jtm_inner_main_exchange" with
   "jgi_jtm_inner_result_queue" tag
6. jtm listens to "jgi_jtm_inner_result_queue" queue; when a result is ready, takes and updates
   tables


Revisions

    05.06.2015 Debugging worker socket disconnection issue;
    05.19.2015 Added channel.close();
    05.20.2015 Added RMQ reconnection for when the connection is lost due to long task processing;

    05.20.2015 2.5.0 released!

    12.03.2015 2.7.0: Tested with heartbeat_interval=0 (RmqConnectionHB(0)) for worker and
                      heartbeat_interval=60 for client in BlockingConnection()

    03.03.2017 3.0.0: Updated to set hearbeat_internal=0 so that connection timeout is disabled
    03.03.2017 3.0.3: Added custom log file location option for worker
    07.10.2017 3.0.4: Updated to print task id when completed

    09.12.2018 3.1.0: Branched out "JTM"; Updated for pika=0.12.0; Changed "type" to "exchange_type"
                      in exchange_declare();
                      Added "jgi_jtm_task_manager" exchange for task and result messages
    09.13.2018 3.1.1: Added unique worker id

    09.13.2018 3.2.0: Added worker heartbeat exchange

    ---------------------------------------------------------------------------------------------

    09.19.2018 0.0.9: Started jtm-worker

    09.21.2018 1.0.0: Working version
    09.27.2018 1.1.0: Updated message structure;

    10.03.2018 1.2.0: Added task_type
    10.04.2018 1.2.1: Set jgi_jtm_client_hb_exchange durable=True and auto_delete=False so that it
                      can be maintained even with no worker
    10.05.2018 1.2.2: Updated resource msg as dict
    10.05.2018 1.2.3: Updated sned hb to client interval to 5sec; Added comm pipe to send taskid
                      with hb;

    10.12.2018 1.4.0: Added recv_reproduce_or_die_proc; jtm_kill done
    10.17.2018 1.4.2: Fixed panic() to terminate child processes
    10.18.2018 1.4.3: Added darwin support for resource reporting
    10.19.2018 1.4.4: Worker can use different queue name (=pool) when user task json has "pool" key
    10.22.2018 1.4.5: Updated workers table; workerId2 for workers table; life_left;
    10.23.2018 1.4.6: Bug fix about workerId2
    10.26.2018 1.4.8: Fixed user termination error code update (-4); Set USER_PROC_PROC_ID = -9 if user
                      process is killed

    10.26.2018 1.5.0: Demo version with static workers tested
    10.28.2018 1.5.1: Added ipaddress to worker hb

    11.06.2018 1.6.0: Tested static workers and sbatch
               1.6.1: Remove user account name from queue name for EC2
    11.08.2018 1.6.2: Added custom queue name postfix for testing in Config
    11.13.2018 1.6.3: Updated jtm_submit for large node sbatch; Updated send_hb to use CNAME as
                      postfix; Added jtm_check_manager cli;

    11.13.2018 1.7.0: Updated sbatch for static worker cloning; removed 'interval' from hb;
    11.14.2018 1.7.1: Added dynamic worker spawning

    11.15.2018 1.8.0: Dynamic workers; Changed hb header format;
    11.19.2018 1.8.1: Changed basic_consume callback args; Changed process_task_request to use the
                      custom pool name as task queue if -p is used; Updated routine to get the
                      current live workers and add a feature to get the num workers per pool name;

    01.16.2019 3.0.1: Fixed options for -wt and -cl;

    02.19.2019 3.2.0: Added Lawrencium slurm support;
    02.25.2019 3.2.3: Set default worker queue name = "small"
    03.11.2019 3.2.8: Fixed invalid child process id in get_pid_tree();
    03.14.2019 3.2.9: Added "--exclusive" for cori sbatch (genepool_special needs it)
    03.21.2019 3.2.10: Added Cori KNL support (KNL constraint doesn't need --qos and --account);

    04.03.2019 3.4.2: Fixed child process killing;

    04.04.2019 4.0.0: Added thread to recv task kill request; Added sleep 10 in sbatch script;

    04.16.2019 4.2.1: Added ack and nack in on_poison(); Added ack and nack in on_kill();
    04.18.2019 4.2.3: Bug fix worker cloning -> basic_reject + requeue=True;

    05.01.2019 5.0.5: Removed scontrol;

    05.09.2019 5.1.1: Fixed life left calculation; Removed raise;

    06.11.2019 5.4.0: Run run_something() with threading.Thread; Added ch._connection.sleep() to fix
                      lost connection issue;
    06.13.2019 5.4.2: Multiprocessing -> threading;
    06.18.2019 5.4.3: Send endData as today in hb to the manager even if the worker is dynamic;

    10.08.2019 5.5.0: Revised recv_hb_from_client_thread;
    10.21.2019 5.5.1: Tested jgi cloud; Updated config; Updated sbatch params;
    10.26.2019 5.6.1: Updated unique worker id + for_loop_index for jaws_lbl_gov;

    02.26.2020 5.6.6: Set slurm job name with pool name;



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

from jaws_jtm.common import setup_custom_logger, logger
from jaws_jtm.config import JtmConfig
from jaws_jtm.lib.rabbitmqconnection import RmqConnectionHB
from jaws_jtm.lib.resourceusage import get_cpu_load, \
    get_runtime, get_pid_tree, get_virtual_memory_usage, \
    get_resident_memory_usage, get_total_mem_usage_per_node, \
    get_num_workers_on_node
from jaws_jtm.lib.run import make_dir, run_sh_command
from jaws_jtm.lib.msgcompress import zdumps, zloads

# This ipc pipe is to send task_id (by on_request())to send_hb_to_client_proc()
# when a task is requested
PIPE_TASK_ID_SEND, PIPE_TASK_ID_RECV = mp.Pipe()

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

config = JtmConfig()
WORKER_TYPE = config.constants.WORKER_TYPE
HB_MSG = config.constants.HB_MSG
VERSION = config.constants.VERSION
COMPUTE_RESOURCES = config.constants.COMPUTE_RESOURCES
TASK_TYPE = config.constants.TASK_TYPE
DONE_FLAGS = config.constants.DONE_FLAGS
NUM_WORKER_PROCS = config.constants.NUM_WORKER_PROCS
TASK_KILL_TIMEOUT_MINUTE = config.constants.TASK_KILL_TIMEOUT_MINUTE

CNAME = config.configparser.get("SITE", "instance_name")
JTM_HOST_NAME = config.configparser.get("SITE", "jtm_host_name")
JTM_INNER_REQUEST_Q = config.configparser.get("JTM", "jtm_inner_request_q")
CTR = config.configparser.getfloat("JTM", "clone_time_rate")
JTM_INNER_MAIN_EXCH = config.configparser.get("JTM", "jtm_inner_main_exch")
JTM_CLIENT_HB_EXCH = config.configparser.get("JTM", "jtm_client_hb_exch")
JTM_WORKER_HB_EXCH = config.configparser.get("JTM", "jtm_worker_hb_exch")
CLIENT_HB_Q_POSTFIX = config.configparser.get("JTM", "client_hb_q_postfix")
WORKER_HB_Q_POSTFIX = config.configparser.get("JTM", "worker_hb_q_postfix")
JTM_TASK_KILL_EXCH = config.configparser.get("JTM", "jtm_task_kill_exch")
JTM_TASK_KILL_Q = config.configparser.get("JTM", "jtm_task_kill_q")
JTM_WORKER_POISON_EXCH = config.configparser.get("JTM", "jtm_worker_poison_exch")
JTM_WORKER_POISON_Q = config.configparser.get("JTM", "jtm_worker_poison_q")
NUM_PROCS_CHECK_INTERVAL = config.configparser.getfloat("JTM", "num_procs_check_interval")
ENV_ACTIVATION = config.configparser.get("JTM", "env_activation")


# -------------------------------------------------------------------------------
def run_something(msg_unzipped, return_msg):
    """
    Run a user command in msg_zipped_to_send
    :param msg_unzipped: uncompressed msg from client
    :param return_msg: msg to return
    :return:
    """
    # Uncompress msg to get a task
    task_id = msg_unzipped["task_id"]
    user_task_cmd = msg_unzipped["user_cmd"]
    out_files = msg_unzipped["output_files"]
    task_type = msg_unzipped["task_type"]
    # cromwell_job_id = msg_unzipped["cromwell_jid"]

    # return_msg = {}
    return_msg["task_id"] = task_id
    return_msg["user_cmd"] = user_task_cmd
    return_msg["task_type"] = task_type

    if "stdout" in msg_unzipped:
        user_task_cmd + " > %s" % (msg_unzipped["stdout"])
    if "stderr" in msg_unzipped:
        user_task_cmd + " 2>%s" % (msg_unzipped["stderr"])

    # Run the task
    logger.info("Running task ID %d...", task_id)
    # logger.info("Running task/task ID, %s/%d...", cromwell_job_id, task_id)

    p = None
    time_out_in_minute = None
    if THIS_WORKER_TYPE != "manual":
        # ex) WORKER_LIFE_LEFT_IN_MINUTE = 20min and TASK_KILL_TIMEOUT_MINUTE = 10min
        # timeout will be set as 10min
        # TASK_KILL_TIMEOUT_MINUTE is a extra housekeeping time after explicitly
        # terminate a task.
        time_out_in_minute = int(WORKER_LIFE_LEFT_IN_MINUTE.value - TASK_KILL_TIMEOUT_MINUTE)
        logger.debug("timeout in minute: %d", time_out_in_minute)

    try:
        p = subprocess.Popen(user_task_cmd.split(),
                             # shell=True,
                             env=os.environ,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
    except Exception as detail:
        logger.exception("Exception: Failed to run user command, %s", user_task_cmd)
        logger.exception("Detail: %s", str(detail))
        raise

    # Set USER_PROC_PROC_ID = forked child process id (this value will be sent
    # to send_hb_to_client function)
    USER_PROC_PROC_ID.value = p.pid
    logger.debug("USER_PROC_PROC_ID.value %d" % USER_PROC_PROC_ID.value)
    logger.debug("User command: %s", user_task_cmd)

    proc_return_code = -1
    stdout_str = None
    time_out_in_second = None
    if time_out_in_minute is not None:
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
            proc_return_code = -6
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
        return_msg["done_flag"] = DONE_FLAGS["failed with user termination"]
        return_msg["ret_msg"] = "Task cancelled."
    else:
        if proc_return_code == -6:
            logger.info("User task timeout. Not enough time left for the worker.")
            return_msg["done_flag"] = "-6"
            return_msg["ret_msg"] = "User task timeout"
        elif proc_return_code == 2:
            logger.info("input file not found.")
            return_msg["done_flag"] = "1"
            return_msg["ret_msg"] = "Input file not found."

        elif proc_return_code == 0:
            logger.info("Task# %s completed!" % (str(task_id)))

            # Output file checking
            if len(out_files):
                ofs = out_files.split(",")
                logger.debug("Number of output files = %d.", len(ofs))
                out_file_list = []

                for i in range(len(ofs)):
                    out_file_list.append(ofs[i])

                FILE_CHECK_INTERVAL = config.configparser.getfloat("JTM", "file_check_interval")
                FILE_CHECKING_MAX_TRIAL = config.configparser.getint("JTM", "file_checking_max_trial")
                FILE_CHECK_INT_INC = config.configparser.getfloat("JTM", "file_check_int_inc")
                ret, file_size = check_output(out_file_list,
                                              FILE_CHECK_INTERVAL,
                                              FILE_CHECKING_MAX_TRIAL,
                                              FILE_CHECK_INT_INC)

                if not ret:
                    ret_msg_str = "Failed to check output file(s): %s, file size = %s." % (ofs, file_size)
                    logger.critical(ret_msg_str)
                    return_msg["done_flag"] = DONE_FLAGS["failed to check output file(s)"]
                    return_msg["ret_msg"] = ret_msg_str
                else:
                    return_msg["done_flag"] = DONE_FLAGS["success with correct output file(s)"]
                    return_msg["ret_msg"] = "Output file checking is OK."
            else:
                return_msg["done_flag"] = DONE_FLAGS["success"]
                return_msg["ret_msg"] = "No file(s) to check."
        else:
            logger.critical("Failed to execute a task, %s. Non-zero exit code. stdout = %s."
                            % (user_task_cmd, stdout_str))
            return_msg["done_flag"] = DONE_FLAGS["failed to run user command"]
            return_msg["ret_msg"] = stdout_str

    # json_data = json.dumps(return_msg)
    logger.debug("Reply msg with result: %s" % str(return_msg))
    # msgZippedToSend = zdumps(json_data)
    # return msgZippedToSend


# -------------------------------------------------------------------------------
def ack_message(ch, delivery_tag):
    """
    :param ch:
    :param delivery_tag:
    :return:
    """
    # Note that `ch` must be the same pika channel instance via which
    #  the message being ACKed was retrieved (AMQP protocol constraint).
    if ch.is_open:
        ch.basic_ack(delivery_tag)
    else:
        # Channel is already closed, so we can't ACK this message;
        # log and/or do something that makes sense for your app in this case.
        pass


# -------------------------------------------------------------------------------
def on_request(ch, method, props, body):
    """
    A callback function whenever a message is received.
    :param ch: channel
    :param method: methodFrame
    :param props: pika connection property
    :param body: message received
    :return:
    """
    msg_unzipped = json.loads(zloads(body))
    task_id = msg_unzipped["task_id"]

    # Get the length of the queue
    # hb_queue_len = ch.queue_declare(queue=method.routing_key,
    #                                 durable=True,
    #                                 exclusive=False,
    #                                 auto_delete=True).method.message_count
    # logger.debug("#tasks queued = %d", hb_queue_len)

    # Send taskid to send_hb_to_client_proc() so that the task_id can be shown in the hb messages
    PIPE_TASK_ID_SEND.send(task_id)

    # If client sends "terminate worker"
    if task_id == -9:
        logger.info("Received TERM signal. Terminate myself.")
        #
        # Ref) http://gavinroy.com/deeper-down-the-rabbit-hole-of-message-redeli
        # Return back the message (the TERM message) in the queue
        # so that the other workers can get it.
        #
        ch.basic_reject(delivery_tag=method.delivery_tag, requeue=True)
        # ch.stop_consuming()
        # ch.close()
        #
        # if recv_hb_from_client_proc_hdl:
        #     recv_hb_from_client_proc_hdl.terminate()
        # if send_hb_to_client_proc_hdl:
        #     send_hb_to_client_proc_hdl.terminate()
        # if recv_poison_proc_hdl:
        #     recv_poison_proc_hdl.terminate()
        # if task_kill_proc_hdl:
        #     task_kill_proc_hdl.terminate()
        #
        # sys.exit(0)
        raise OSError(2, 'Worker termination request')

    logger.info("Received a task, %r" % (msg_unzipped,))
    logger.debug("Return queue = %s", props.reply_to)

    ######################################
    # OLD
    # response = run_something(msg_unzipped)

    # NEW
    # Note: to fix connection lost
    # ex) 2019-06-10 12:13:37,917 | jtm-worker | on_request | CRITICAL : Something wrong in on_request():
    # Stream connection lost: error(104, 'Connection reset by peer')
    #
    # https://stackoverflow.com/questions/14572020/handling-long-running-tasks-in-pika-rabbitmq
    #
    result_dict = {}
    thread = threading.Thread(target=run_something,
                              args=(msg_unzipped, result_dict))
    thread.start()
    while thread.is_alive():  # Loop while the thread is processing
        ch._connection.sleep(1.0)
    json_data = json.dumps(result_dict)
    logger.debug("Reply msg with result: %s" % str(json_data))
    logger.debug(json_data)
    response = zdumps(json_data)
    ######################################

    try:
        logger.debug("Sent the result back to the client via '%s' queue.", props.reply_to)

        # Note: to keep result messages even the client is not alive.
        ch.queue_declare(queue=props.reply_to,
                         durable=True,
                         exclusive=False,
                         auto_delete=True)

        assert props.reply_to.endswith(CNAME)
        ch.basic_publish(exchange=JTM_INNER_MAIN_EXCH,
                         routing_key=props.reply_to,  # use the queue which the client created
                         properties=pika.BasicProperties(
                             delivery_mode=2,  # make message persistent
                             correlation_id=props.correlation_id),
                         body=response)

        # Note: After sending ack, the message will be deleted from RabbitMQ
        #  If this worker crashes while running a user command, this task will
        #  be sent to other workers available
        ch.basic_ack(delivery_tag=method.delivery_tag)

    except Exception as e:
        logger.critical("Something wrong in on_request(): %s", e)
        raise

    # Send taskid=0 to send_hb_to_client_proc() b/c the requested task is completed
    PIPE_TASK_ID_SEND.send(0)

    # Reset child pid
    USER_PROC_PROC_ID.value = 0


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

                if b_is_file_exist and file_size > 0:
                    b_is_file_found = True
                    logger.info("Output file '%s' is OK.", a_file)
                    break
            else:
                logger.info("Outout file not found.")

            # Wait for initial wait time
            time.sleep(out_file_check_wait_time)
            # conn.process_data_events(time_limit=float(worker_timeout_in_sec))

            # Increase the wait time
            out_file_check_wait_time *= out_file_check_wait_time_increase
            trial += 1

    return b_is_file_found, file_size


# -------------------------------------------------------------------------------
def send_hb_to_client_proc(interval, pipe_task_id,
                           slurm_job_id, mem_per_node, mem_per_core, num_cores,
                           job_time, clone_time_rate, task_queue_name, pool_name, nwpn,
                           exch_name, worker_hb_queue):
    """
    Send heartbeats to the client
    :param interval: time interval to send heartbeats to the client
    :param pipe_task_id: Pythin IPC pipe for getting task id
    :param slurm_job_id: SLURM job id
    :param mem_per_node: memory request per node
    :param mem_per_core: memory request per core
    :param num_cores: number of cores
    :param job_time: wallclocktime
    :param clone_time_rate: clone time rate
    :param task_queue_name: task queue name
    :param pool_name: pool name
    :param nwpn: number of workers per node
    :return:
    """
    # Remote broker (mq.nersc.gov)
    # rmq_conn = RmqConnection()
    # with heartbeat_interval=0
    # ref) http://stackoverflow.com/questions/14572020/handling-long-running-tasks-in-pika-rabbitmq,
    # http://stackoverflow.com/questions/34721178/pika-blockingconnection-rabbitmq-connection-closed
    rmq_conn = RmqConnectionHB()
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

    while 1:
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
            #
            logger.debug("rootpid = {} childpid = {}".format(PARENT_PROCESS_ID, USER_PROC_PROC_ID.value))
            if USER_PROC_PROC_ID.value == 0:
                child_pid = PARENT_PROCESS_ID
            else:
                child_pid = int(USER_PROC_PROC_ID.value)

            # Collect pids from process tree
            root_pid_num = get_pid_tree(PARENT_PROCESS_ID)
            logger.debug("get_pid_tree(root_proc_id={}) = {}".format(PARENT_PROCESS_ID, root_pid_num))

            if PARENT_PROCESS_ID != child_pid:
                if child_pid == -9:  # if it's terminated by "kill"
                    child_pid = PARENT_PROCESS_ID
                pid_list_child = get_pid_tree(child_pid)
                logger.debug("get_pid_tree(child_pid={}) = {}".format(child_pid, pid_list_child))
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
                cpu_load_list = [float(get_cpu_load(pid)) for pid in proc_id_list_merged]
            except Exception as e:
                logger.exception("get_cpu_load() exception: {}".format(e))
                # ch.stop_consuming()
                # ch.close()
                # conn.close()
                # os._exit(1)
                raise

            max_cpu_load = max(cpu_load_list) if len(cpu_load_list) > 0 else 0.0

            # Collect mem_usages for all pids in the tree and get sum()
            rmem_usage = "%.1f" % sum(rmem_usage_list)
            vmem_usage = "%.1f" % sum(vmem_usage_list)

            # TODO
            # Compare mem_per_node with rmem_usage or vmem_usage
            # if mem usage >= mem_per_node, terminate the task (=: margin)
            # using

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
            num_workers_on_node = get_num_workers_on_node()

            # Check if there is any task id in the ipc pipe
            if pipe_task_id.poll():
                task_id = pipe_task_id.recv()

            today = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            end_date_time = None
            # worker_life_left_in_minute = 0
            # WORKER_LIFE_LEFT_IN_MINUTE.value = 0

            try:
                # Note: this only works with SLURM. For other HPCs/Clouds, methods to get the
                #  remaining life are needed! ==> resolved!
                if slurm_job_id:
                    # OLD
                    # fixme; slurm overload
                    # so, se, ec = run_sh_command("scontrol show jobid %d"
                    #                             % (slurm_job_id), live=True, log=logger, stdoutPrint=False)
                    # pat = re.compile(r"""([^\s=]+)=\s*((?:[^\s=]+(?:\s|$))*)""")
                    # # Convert a list of x=y to dict
                    # entries = dict((k, v.split()) for k, v in pat.findall(so))
                    # delta = datetime.datetime.now() - datetime.datetime.strptime(entries["EndTime"][0],
                    # "%Y-%m-%dT%H:%M:%S")
                    # worker_life_left_in_minute = -divmod(delta.total_seconds(), 60)[0]
                    # end_date_time = entries["EndTime"][0]

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
                    #  = divmod(delta.total_seconds(), 60)[0]
                    WORKER_LIFE_LEFT_IN_MINUTE.value = int(divmod(delta.total_seconds(), 60)[0])
            except Exception as e:
                logger.critical("Something wrong in computing remaining wall clock time: %s", e)
                # ch.stop_consuming()
                # ch.close()
                # conn.close()
                # os._exit(1)
                raise

            # Get the total number of tasks in the task queue for this pool
            hb_queue_len = ch.queue_declare(queue=task_queue_name,
                                            durable=True,
                                            exclusive=False,
                                            auto_delete=True).method.message_count
            # logger.debug("#tasks queued=%d in %s" % (hb_queue_len, task_queue_name))

            # To decrease the hb message size, changed to int key value
            try:
                ip_address = socket.gethostbyname(host_name)
            except Exception:
                # eprint("Exception: Failed to get ip address: %s" % e)
                ip_address = None

            ncore_param = mp.cpu_count()
            if WORKER_TYPE[THIS_WORKER_TYPE] > 0:
                ncore_param = num_cores

            msg_dict_to_send = {HB_MSG["child_pid"]: child_pid,
                                HB_MSG["clone_time_rate"]: clone_time_rate,
                                HB_MSG["cpu_load"]: max_cpu_load,
                                HB_MSG["end_date"]: today,  # Note: for dynamic worker endDate update
                                HB_MSG["host_name"]: host_name,
                                HB_MSG["ip_address"]: ip_address,
                                HB_MSG["job_time"]: job_time if WORKER_TYPE[THIS_WORKER_TYPE] > 0 else None,
                                HB_MSG["jtm_host_name"]: JTM_HOST_NAME,
                                HB_MSG["life_left"]: WORKER_LIFE_LEFT_IN_MINUTE.value,
                                HB_MSG["mem_per_core"]: mem_per_core if WORKER_TYPE[THIS_WORKER_TYPE] > 0 else "",
                                HB_MSG["mem_per_node"]: mem_per_node if WORKER_TYPE[THIS_WORKER_TYPE] > 0 else "",
                                HB_MSG["num_cores"]: ncore_param,
                                HB_MSG["num_tasks"]: hb_queue_len,
                                HB_MSG["num_workers_on_node"]: num_workers_on_node,
                                HB_MSG["perc_mem_used"]: perc_used_mem,
                                HB_MSG["pool_name"]: pool_name,
                                HB_MSG["ret_msg"]: "hb",
                                HB_MSG["rmem_usage"]: rmem_usage,
                                HB_MSG["root_pid"]: PARENT_PROCESS_ID,
                                HB_MSG["run_time"]: proc_run_time,
                                HB_MSG["slurm_jobid"]: slurm_job_id,
                                HB_MSG["task_id"]: task_id,
                                HB_MSG["vmem_usage"]: vmem_usage,
                                HB_MSG["worker_id"]: UNIQ_WORKER_ID,
                                HB_MSG["worker_type"]: WORKER_TYPE[THIS_WORKER_TYPE],
                                HB_MSG["nwpn"]: nwpn}

            msg_zipped_to_send = zdumps(json.dumps(msg_dict_to_send))

            assert worker_hb_queue.endswith(CNAME)
            ch.basic_publish(exchange=exch_name,
                             routing_key=worker_hb_queue,
                             body=msg_zipped_to_send)

            logger.debug("Send HB to the client: {}".format(msg_dict_to_send))

        except Exception as e:
            logger.critical("Something wrong with send_hb_to_client(): %s", e)
            # Note: Need to terminate the whole program, not this thread only.
            # ch.stop_consuming()
            # ch.close()
            # conn.close()
            # os._exit(1)
            raise

        # Todo: Need to start with shorted interval like (0.5-1 sec) for about 30 sec
        #  from the beginning and then use the user specified interval
        # time.sleep(float(interval))
        conn.process_data_events(time_limit=float(interval))

    # unreachable
    ch.close()
    conn.close()


# -----------------------------------------------------------------------
def recv_hb_from_client_proc2(task_queue, exch_name, cl_hb_q_postfix):
    """

    :param task_queue:
    :param exch_name:
    :param cl_hb_q_postfix:
    :return:
    """
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()
    # host_name = socket.gethostname()

    # Declare exchange
    ch.exchange_declare(exchange=exch_name,
                        # exchange_type="fanout",
                        exchange_type="topic",
                        durable=False,
                        auto_delete=False)

    # Declare queue
    client_hb_queue_name = "_jtm_worker_%s%s" % (UNIQ_WORKER_ID, cl_hb_q_postfix)
    ch.queue_declare(queue=client_hb_queue_name,
                     durable=False,
                     exclusive=True,
                     auto_delete=True)

    # Bind exchange to hb queue for receiving the client hb
    ch.queue_bind(exchange=exch_name,
                  queue=client_hb_queue_name,
                  routing_key="*." + CNAME)

    logger.info('[*] Waiting for the manager.')

    # def stop():
    #     logger.info("Terminate myself: {} on {}".format(UNIQ_WORKER_ID, host_name))
    #     if recv_hb_from_client_proc_hdl:
    #         recv_hb_from_client_proc_hdl.terminate()
    #     if send_hb_to_client_proc_hdl:
    #         send_hb_to_client_proc_hdl.terminate()
    #     if recv_poison_proc_hdl:
    #         recv_poison_proc_hdl.terminate()
    #     if task_kill_proc_hdl:
    #         task_kill_proc_hdl.terminate()
    #     ch.stop_consuming()
    #     # os.kill(int(ppid), signal.SIGTERM)  # term jtm-worker process
    #     # sys.exit(10)
    #     ch.close()
    #     conn.close()
    #     os._exit(1)

    def callback(ch, method, properties, body):
        msg_unzipped = json.loads(zloads(body))
        if msg_unzipped["task_type"] == TASK_TYPE["term"] and msg_unzipped["task_queue"] == task_queue:
            # stop()
            raise OSError(2, 'Worker terminatin request received.')
        elif msg_unzipped["task_type"] == TASK_TYPE["hb"]:
            # Poison from the manager. Terminate itself.
            logger.debug("Received heartbeats from the client, %r:%r" % (method.routing_key, msg_unzipped))

    ch.basic_consume(queue=client_hb_queue_name,
                     on_message_callback=callback,
                     auto_ack=True)

    ch.start_consuming()


# -------------------------------------------------------------------------------
def recv_task_kill_request_proc():
    """
    Wait for task termination request from JTM

    :param exch_name:
    :param queue_name:
    :return:
    """

    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()
    exch_name = JTM_TASK_KILL_EXCH
    queue_name = JTM_TASK_KILL_Q
    routing_key = UNIQ_WORKER_ID

    ch.exchange_declare(exchange=exch_name,
                        exchange_type="fanout",
                        durable=True,
                        auto_delete=False)
    ch.queue_declare(queue=queue_name,
                     # durable=False,
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

                # This -9 is to notify run_something() that it's killed by user requests
                # Also send_hb_to_client() will check this for adjust childpid and parentpid
                # Note: this should done first to signal run_something() that the process is killed.
                USER_PROC_PROC_ID.value = -9

                # kill if there is child's children
                for i in get_pid_tree(msg_unzipped["child_pid"]):
                    # kill_cmd = "kill -9 %d" % msg_unzipped["child_pid"]
                    kill_cmd = "kill -9 %d" % i
                    logger.info("Executing {} for taskID, {}".format(kill_cmd, msg_unzipped["task_id"]))
                    _, _, ec = run_sh_command(kill_cmd, live=True, log=logger)
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
                _, _, ec = run_sh_command(kill_cmd, live=True, log=logger)
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
def recv_reproduce_or_die_proc(task_queue, cluster_name, mem_per_node, mem_per_core,
                               num_nodes, num_cores, wallclocktime, clone_time_rate,
                               nwpn, exch_name, queue_name):
    """
    Wait for process termination request from JTM
    :param task_queue: task queue (pool) name
    :param cluster_name: cluster name
    :param mem_per_node: memory size per node
    :param mem_per_core: memory size per core
    :param num_nodes: number of nodes
    :param num_cores: number of cores
    :param wallclocktime: wallclocktime request
    :param clone_time_rate: cloning time rate (remaining runtime / total runtime)
    :param nwpn: number of workers per node
    :return:
    """
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()

    ch.exchange_declare(exchange=exch_name,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)
    ch.queue_declare(queue=queue_name,
                     # durable=False,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)
    ch.queue_bind(exchange=exch_name,
                  queue=queue_name,
                  routing_key=UNIQ_WORKER_ID)

    def on_poison(ch, method, props, body):
        msg_unzipped = json.loads(zloads(body))

        if msg_unzipped["worker_id"] == UNIQ_WORKER_ID:

            # If static worker cloning request
            if msg_unzipped["num_clones"] == -1:
                # Kill the worker request
                pass
            elif msg_unzipped["num_clones"] == 1:
                logger.info("Static worker cloning request received.")

                jtm_worker_cmd = "jtm worker -wt static -cl {} -t {} -ct {}".format(
                                  cluster_name, wallclocktime, clone_time_rate)

                if task_queue is not None:
                    jtm_worker_cmd += " -p {}".format(task_queue)

                if nwpn != "" and int(nwpn) > 1:
                    jtm_worker_cmd += " -nw {}".format(nwpn)
                else:
                    jtm_worker_cmd += " -nw 1"

                if num_nodes != "" and int(num_nodes) > 0:
                    jtm_worker_cmd += " -N {} -m {} -c {}".format(num_nodes, mem_per_node, num_cores)
                else:
                    jtm_worker_cmd += " -c {}".format(num_cores)
                    if mem_per_node not in ("", None):  # if node mem defined, set --mem
                        jtm_worker_cmd += " -m {}".format(mem_per_node)
                    else:  # if not, set mem per core
                        jtm_worker_cmd += " -mc {}".format(mem_per_core)

                logger.info("Executing {}".format(jtm_worker_cmd))
                so, se, ec = run_sh_command(jtm_worker_cmd, live=True, log=logger)
                logger.debug("{} {} {}".format(so, se, ec))
                run_sh_command(jtm_worker_cmd + " --dry-run", live=True, log=logger)

            ch.basic_ack(delivery_tag=method.delivery_tag)

        else:
            # ch.basic_nack(delivery_tag=method.delivery_tag)
            # If UNIQ_WORKER_ID is not for me, reject and requeue it
            ch.basic_reject(delivery_tag=method.delivery_tag,
                            requeue=True)

    # waiting for poison or cloning command
    ch.basic_qos(prefetch_count=1)
    ch.basic_consume(queue=queue_name,
                     on_message_callback=on_poison,
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

    """
    # ps_cmd = "ps -aef | grep '%s' | grep %s | grep -v grep | awk '{print $2}'" \
    #          % ("jtm worker", getpass.getuser())
    # if sys.platform.lower() == "darwin":
    #     ps_cmd = "ps -aef | grep '%s' | grep -v grep | awk '{print $2}'" \
    #              % ("jtm worker")
    # while 1:
    #     pid_list = back_ticks(ps_cmd, shell=True)
    #     if type(pid_list) is bytes:
    #         pid_list = pid_list.decode()
    #     pid_list = [int(i) for i in pid_list.split('\n') if i]
    #     logger.debug("ps_cmd: {}".format(ps_cmd))
    #     logger.debug("pid list = {}".format(pid_list))
    #     logger.debug("len(pid) = {}".format(len(pid_list)))
    #     # if len(pid_list) != NUM_WORKER_PROCS:
    #     #     logger.critical("JTM worker child process error.")
    #     #     [os.kill(i, signal.SIGTERM) for i in pid_list]
    #     time.sleep(NUM_PROCS_CHECK_INTERVAL)

    while 1:
        logger.debug("Proc ID list: {}".format(pid_list))
        for pid in pid_list:
            try:
                os.kill(pid, 0)
            except OSError:
                logger.critical("Child process doesn't exist: %d" % (pid))
                raise

        time.sleep(NUM_PROCS_CHECK_INTERVAL)


# -------------------------------------------------------------------------------
def worker(heartbeat_interval_param, custom_log_dir, custom_job_log_dir_name,
           pool_name_param, worker_timeout_in_sec_param: int, dry_run: bool,
           slurm_job_id_param: int, worker_type_param, cluster_name_param,
           worker_clone_time_rate_param: float, num_workers_per_node_param: int,
           worker_id_param: str, charging_account_param, num_nodes_to_request_param: int,
           num_cores_to_request_param: int, constraint_param,
           mem_per_node_to_request_param: str, mem_per_cpu_to_request_param: str,
           qos_param, job_time_to_request_param: str, debug: bool) -> int:

    # Job dir setting
    job_script_dir_name = os.path.join(config.configparser.get("JTM", "log_dir"), "job")
    if custom_job_log_dir_name:
        job_script_dir_name = custom_job_log_dir_name
    make_dir(job_script_dir_name)

    # Log dir setting
    log_dir_name = os.path.join(config.configparser.get("JTM", "log_dir"), "log")
    if custom_log_dir:
        log_dir_name = custom_log_dir
    make_dir(log_dir_name)

    print("JTM Worker, version: {}".format(VERSION))

    # Set uniq worker id if worker id is provided in the params
    if worker_id_param:
        global UNIQ_WORKER_ID
        UNIQ_WORKER_ID = worker_id_param

    # Logger setting
    log_level = "info"
    if debug:
        log_level = "debug"

    setup_custom_logger(log_level, log_dir_name,
                        1, 1,
                        worker_id=UNIQ_WORKER_ID)

    logger.info("\n*****************\nDebug mode is %s\n*****************"
                % ("ON" if debug else "OFF"))

    RMQ_HOST = config.configparser.get("RMQ", "host")
    RMQ_PORT = config.configparser.get("RMQ", "port")
    USER_NAME = config.configparser.get("SITE", "user_name")
    PRODUCTION = False
    if config.configparser.get("JTM", "run_mode") == "prod":
        PRODUCTION = True
    JOBTIME = config.configparser.get("SLURM", "jobtime")
    CONSTRAINT = config.configparser.get("SLURM", "constraint")
    CHARGE_ACCNT = config.configparser.get("SLURM", "charge_accnt")
    QOS = config.configparser.get("SLURM", "qos")
    PARTITION = config.configparser.get("SLURM", "partition")
    MEMPERCPU = config.configparser.get("SLURM", "mempercpu")
    NWORKERS = config.configparser.getint("JTM", "num_workers_per_node")
    NCPUS = config.configparser.getint("SLURM", "ncpus")

    # # Todo: site specific setting --> remove
    # CORI_KNL_CHARGE_ACCNT = config.configparser.get("SLURM", "knl_charge_accnt")
    # CORI_KNL_QOS = config.configparser.get("SLURM", "knl_qos")
    hearbeat_interval = config.configparser.getfloat("JTM", "worker_hb_send_interval")
    worker_timeout_in_sec = config.configparser.getfloat("JTM", "worker_timeout")

    logger.info("Set jtm log file location to %s", log_dir_name)
    logger.info("Set jtm job file location to %s", job_script_dir_name)
    logger.info("RabbitMQ broker: %s", RMQ_HOST)
    logger.info("RabbitMQ port: %s", RMQ_PORT)
    logger.info("Pika version: %s", pika.__version__)
    logger.info("JTM user name: %s", USER_NAME)
    logger.info("Unique worker ID: %s", UNIQ_WORKER_ID)
    logger.info("\n*****************\nRun mode is %s\n*****************"
                % ("PROD" if PRODUCTION else "DEV"))
    logger.debug("env activation: %s", ENV_ACTIVATION)

    # Slurm config
    num_nodes_to_request = 0
    if num_nodes_to_request_param:
        num_nodes_to_request = num_nodes_to_request_param
        # Todo
        #  Cori and JGI Cloud are exclusive allocation. So this is not needed.
        # assert mem_per_node_to_request_param is not None, "-N needs --mem-per-cpu (-mc) setting."

    ###########################################################################
    # 11.13.2018 decided to remove all default values from argparse
    num_workers_per_node = num_workers_per_node_param if num_workers_per_node_param else NWORKERS
    num_cpus_to_request = num_cores_to_request_param if num_cores_to_request_param else NCPUS
    mem_per_cpu_to_request = mem_per_cpu_to_request_param if mem_per_cpu_to_request_param else MEMPERCPU
    mem_per_node_to_request = mem_per_node_to_request_param if mem_per_node_to_request_param else ""
    ###########################################################################

    job_time_to_request = job_time_to_request_param if job_time_to_request_param else JOBTIME

    constraint = constraint_param if constraint_param else CONSTRAINT
    charging_account = charging_account_param if charging_account_param else CHARGE_ACCNT
    qos = qos_param if qos_param else QOS

    global THIS_WORKER_TYPE
    THIS_WORKER_TYPE = worker_type_param
    job_name = "jtm_worker_" + pool_name_param

    # Set task queue name
    inner_task_request_queue = None
    if heartbeat_interval_param:
        hearbeat_interval = heartbeat_interval_param
    if worker_timeout_in_sec_param:
        worker_timeout_in_sec = worker_timeout_in_sec_param

    # Start hb receive thread
    tp_name = ""
    if pool_name_param:
        tp_name = pool_name_param

    # If you want to create a custom worker pool, use this.
    # JTM detects "pool" field from task json, and creates custom pool if it's set.
    # All tasks with the pool name will be directed to the custom pool
    # user_account_name = getpass.getuser().lower()

    # Todo: feature to add a manual worker to a specific cluster??
    #  ex) $ jtm-worker -cl cori -p small (on cori)
    # if worker_type_param == "manual" and cluster_name_param != "local":
    #     # inner_task_request_queue = "_jtm_inner_request_queue." + cluster_name_param + "." + pool_name_param
    #     JTM_INNER_REQUEST_Q = "_jtm_inner_request_queue." + cluster_name_param + "." + USER_NAME

    assert pool_name_param is not None, "User pool name is not set"
    inner_task_request_queue = JTM_INNER_REQUEST_Q + "." + pool_name_param
    # else:  # not reachable. default inner_task_request_queue = small
    #     inner_task_request_queue = JTM_INNER_REQUEST_Q

    worker_clone_time_rate = worker_clone_time_rate_param if worker_clone_time_rate_param else CTR
    if THIS_WORKER_TYPE in ("static", "dynamic"):
        assert cluster_name_param != "" and \
               cluster_name_param != "local", "Static or dynamic worker needs a cluster setting (-cl)."

    slurm_job_id = slurm_job_id_param
    cluster_name = cluster_name_param

    if cluster_name == "cori" and mem_per_cpu_to_request != "" and \
            float(mem_per_cpu_to_request.replace("GB", "").replace("G", "").replace("gb", "")) > 1.0:
        logger.critical("--mem-per-cpu in Cori shouldn't be larger than 1GB. User '--mem' instead.")
        sys.exit(1)

    logger.info("RabbitMQ broker: %s", RMQ_HOST)
    logger.info("Task queue name: %s", inner_task_request_queue)
    logger.info("Worker type: %s", THIS_WORKER_TYPE)

    if slurm_job_id == 0 and THIS_WORKER_TYPE in ["static", "dynamic"]:
        batch_job_script_file = os.path.join(job_script_dir_name, "jtm_%s_worker_%s.job" %
                                             (THIS_WORKER_TYPE, UNIQ_WORKER_ID))
        batch_job_script_str = ""
        batch_job_misc_params = ""

        if cluster_name in ("cori", "lawrencium", "jgi_cloud", "jaws_lbl_gov", "lbl", "jgi_cluster"):

            with open(batch_job_script_file, "w") as jf:
                batch_job_script_str += "#!/bin/bash -l"

                if cluster_name in ("cori"):

                    if num_nodes_to_request_param:
                        batch_job_script_str += """
#SBATCH -N %(num_nodes_to_request)d
#SBATCH --mem=%(mem)s""" % dict(num_nodes_to_request=num_nodes_to_request, mem=mem_per_node_to_request)
                        batch_job_misc_params += " -N %(num_nodes_to_request)d -m %(mem)s" % \
                                                 dict(num_nodes_to_request=num_nodes_to_request,
                                                      mem=mem_per_node_to_request)

                        if num_cores_to_request_param:
                            batch_job_script_str += """
#SBATCH -c %(num_cores)d""" % dict(num_cores=num_cpus_to_request)
                            batch_job_misc_params += " -c %(num_cores)d" % \
                                                     dict(num_cores=num_cpus_to_request)

                    else:
                        batch_job_script_str += """
#SBATCH -c %(num_cores)d""" % dict(num_cores=num_cpus_to_request)
                        batch_job_misc_params += " -c %(num_cores)d" % \
                                                 dict(num_cores=num_cpus_to_request)

                        if mem_per_node_to_request_param:
                            batch_job_script_str += """
#SBATCH --mem=%(mem)s""" % dict(mem=mem_per_node_to_request)
                            batch_job_misc_params += " -m %(mem)s " % \
                                                     dict(mem=mem_per_node_to_request)
                        else:
                            batch_job_script_str += """
#SBATCH --mem-per-cpu=%(mempercore)s""" % dict(mempercore=mem_per_cpu_to_request)
                            batch_job_misc_params += " -mc %(mempercore)s" % \
                                                     dict(mempercore=mem_per_cpu_to_request)

                        if worker_id_param:
                            batch_job_misc_params += " -wi %(worker_id)s${i}" % \
                                                     dict(worker_id=UNIQ_WORKER_ID)

                    ###########################
                    if cluster_name == "cori":

                        # Need to set both --qos=genepool (or genepool_shared) _and_ -A fungalp
                        # OR
                        # no qos _and_ -A m342 _and_ -C haswell

                        # Note: currently constraint in ["haswell" | "knl"]
                        if constraint == "haswell":
                            if qos_param:
                                batch_job_script_str += """
#SBATCH -q %(qosname)s""" % dict(qosname=qos)
                                batch_job_misc_params += " -q %(qosname)s" % dict(qosname=qos)

                            else:
                                batch_job_script_str += """
#SBATCH -q %(qosname)s""" % dict(qosname=qos)

                            batch_job_script_str += """
#SBATCH -C haswell"""

                            if charging_account == "m342":
                                batch_job_misc_params += " -A %(sa)s" % dict(sa="m342")

                            batch_job_script_str += """
#SBATCH -A %(charging_account)s""" % dict(charging_account=charging_account)

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
                            batch_job_script_str += """
#SBATCH -C knl
#SBATCH -A %(charging_account)s
#SBATCH -q %(qosname)s""" % \
                                                    dict(charging_account=charging_account, qosname=qos)

                            batch_job_misc_params += " -A %(charging_account)s -q %(qosname)s" % \
                                                     dict(charging_account=charging_account, qosname=qos)

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

                            batch_job_script_str += """
#SBATCH -C skylake
#SBATCH -A %(charging_account)s
#SBATCH -q %(qosname)s""" % \
                                                    dict(charging_account=charging_account, qosname=qos)

                            batch_job_misc_params += " -A %(charging_account)s -q %(qosname)s" % \
                                                     dict(charging_account=charging_account, qosname=qos)

                        excl_param = ""
                        if constraint != "skylake":
                            excl_param = "#SBATCH --exclusive"

                        tq_param = ""
                        if pool_name_param:
                            tq_param = "-p " + pool_name_param

                        batch_job_script_str += """
#SBATCH -t %(wall_time)s
#SBATCH --job-name=%(job_name)s
#SBATCH -o %(job_dir)s/jtm_%(worker_type)s_worker_%(worker_id)s.out
#SBATCH -e %(job_dir)s/jtm_%(worker_type)s_worker_%(worker_id)s.err
%(exclusive)s

module unload python
%(env_activation_cmd)s
for i in {1..%(num_workers_per_node)d}
do
    echo "jobid: $SLURM_JOB_ID"
    jtm worker %(debug)s --slurm_job_id $SLURM_JOB_ID \
-cl cori \
-wt %(worker_type)s \
-t %(wall_time)s \
--clone_time_rate %(clone_time_rate)f %(task_queue)s \
--num_worker_per_node %(num_workers_per_node)d \
-C %(constraint)s \
%(other_params)s &
    sleep 1
done
wait
""" % \
                                                dict(debug="--debug" if debug else "",
                                                     wall_time=job_time_to_request,
                                                     job_dir=job_script_dir_name,
                                                     worker_id=UNIQ_WORKER_ID,
                                                     worker_type=THIS_WORKER_TYPE,
                                                     clone_time_rate=worker_clone_time_rate,
                                                     task_queue=tq_param,
                                                     num_workers_per_node=num_workers_per_node,
                                                     env_activation_cmd=ENV_ACTIVATION,
                                                     other_params=batch_job_misc_params,
                                                     constraint=constraint,
                                                     job_name=job_name,
                                                     exclusive=excl_param)

                elif cluster_name in ("lawrencium", "jgi_cloud", "jaws_lbl_gov", "jgi_cluster", "lbl"):

                    if worker_id_param:
                        # batch_job_misc_params += " -wi %s" % (worker_id_param)
                        batch_job_misc_params += " -wi %(worker_id)s${i}" \
                                                 % dict(worker_id=UNIQ_WORKER_ID)

                    tp_param = ""
                    if pool_name_param:
                        tp_param = "-p " + pool_name_param
                    part_param = ""
                    if cluster_name == "lawrencium":
                        part_param = PARTITION
                    else:
                        part_param = PARTITION
                    qos_param = ""
                    if cluster_name == "lawrencium":
                        qos_param = QOS
                    else:
                        qos_param = QOS
                    charge_param = ""
                    if cluster_name == "lawrencium":
                        charge_param = CHARGE_ACCNT
                    else:
                        charge_param = CHARGE_ACCNT
                    nnode_param = 1
                    if num_nodes_to_request_param:
                        nnode_param = num_nodes_to_request

                    mnode_param = "#SBATCH --mem=%(mem)s" \
                                  % dict(mem=mem_per_node_to_request)

                    setenv_param = ""
                    if cluster_name == "jaws_lbl_gov":
                        setenv_param = "export JTM_HOST_NAME=jaws_lbl_gov"

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
%(set_env_for_jtm_host_name)s
%(env_activation_cmd)s

for i in {1..%(num_workers_per_node)d}
do
    echo "jobid: $SLURM_JOB_ID"
    jtm worker %(debug)s --slurm_job_id $SLURM_JOB_ID \
-cl %(lbl_cluster_name)s \
-wt %(worker_type)s \
-t %(wall_time)s \
--clone_time_rate %(clone_time_rate)f %(task_queue)s \
--num_worker_per_node %(num_workers_per_node)d \
%(other_params)s &
    sleep 1
done
wait
""" % \
                                            dict(debug="--debug" if debug else "",
                                                 wall_time=job_time_to_request,
                                                 job_name=job_name,
                                                 partition_name=part_param,
                                                 qosname=qos_param,
                                                 charging_account=charge_param,
                                                 num_nodes_to_request=nnode_param,
                                                 mem_per_node_setting=mnode_param,
                                                 worker_id=UNIQ_WORKER_ID,
                                                 job_dir=job_script_dir_name,
                                                 set_env_for_jtm_host_name=setenv_param,
                                                 env_activation_cmd=ENV_ACTIVATION,
                                                 num_workers_per_node=num_workers_per_node,
                                                 lbl_cluster_name=cluster_name,
                                                 worker_type=THIS_WORKER_TYPE,
                                                 clone_time_rate=worker_clone_time_rate,
                                                 task_queue=tp_param,
                                                 other_params=batch_job_misc_params)

                jf.writelines(batch_job_script_str)

            os.chmod(batch_job_script_file, 0o775)

            if dry_run:
                print(batch_job_script_str)
                sys.exit(0)

            _, _, ec = run_sh_command("sbatch --parsable %s"
                                      % (batch_job_script_file), live=True, log=logger)
            assert ec == 0, "Failed to run 'jtm worker' to sbatch dynamic worker."
            return ec

        elif cluster_name == "aws":
            pass

    # If it's spawned by sbatch
    # Todo: need to record job_id, worker_id, worker_type, starting_time, wallclocktime
    # scontrol show jobid -dd <jobid> ==> EndTime
    # scontrol show jobid <jobid> ==> EndTime
    # sstat --format=AveCPU,AvePages,AveRSS,AveVMSize,JobID -j <jobid> --allsteps
    #
    # if endtime - starttime <= 10%, execute sbatch again
    # if slurm_job_id != 0 and THIS_WORKER_TYPE == "static":
    #     logger.debug("worker_type: {}".format(THIS_WORKER_TYPE))
    #     logger.debug("slurm_job_id: {}".format(slurm_job_id))

    # Dynamic workers creates [[two]] children when it approaches to the wallclocktime limit
    # considering the task queue length
    # Also, maintain the already requested number of workers
    # if no more workers needed, it won't call sbatch
    # elif slurm_job_id != 0 and THIS_WORKER_TYPE == "dynamic":
    #     logger.debug("worker_type: {}".format(THIS_WORKER_TYPE))
    #     logger.debug("slurm_job_id: {}".format(slurm_job_id))

    # Remote broker (rmq.nersc.gov)
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()
    # ch.confirm_delivery()
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
    logger.debug("Main pid = {}".format(PARENT_PROCESS_ID))

    pid_list = []
    pid_list.append(PARENT_PROCESS_ID)

    def proc_clean():
        if recv_hb_from_client_proc_hdl:
            recv_hb_from_client_proc_hdl.terminate()
        if send_hb_to_client_proc_hdl:
            send_hb_to_client_proc_hdl.terminate()
        if recv_poison_proc_hdl:
            recv_poison_proc_hdl.terminate()
        if task_kill_proc_hdl:
            task_kill_proc_hdl.terminate()
        if check_processes_hdl:
            check_processes_hdl.terminate()

    # Start task termination thread
    task_kill_proc_hdl = mp.Process(target=recv_task_kill_request_proc)
    task_kill_proc_hdl.daemon = True
    task_kill_proc_hdl.start()
    pid_list.append(task_kill_proc_hdl.pid)

    try:
        recv_hb_from_client_proc_hdl = mp.Process(target=send_hb_to_client_proc,
                                                  args=(hearbeat_interval,
                                                        PIPE_TASK_ID_RECV,
                                                        slurm_job_id,
                                                        mem_per_node_to_request,
                                                        mem_per_cpu_to_request,
                                                        num_cpus_to_request,
                                                        job_time_to_request,
                                                        worker_clone_time_rate,
                                                        inner_task_request_queue,
                                                        tp_name,
                                                        num_workers_per_node,
                                                        JTM_WORKER_HB_EXCH,
                                                        WORKER_HB_Q_POSTFIX))

        recv_hb_from_client_proc_hdl.daemon = True
        recv_hb_from_client_proc_hdl.start()
    except Exception as e:
        logger.exception("send_hb_to_client_proc: {}".format(e))
        proc_clean()
        ch.stop_consuming()
        ch.close()
        conn.close()
        sys.exit(1)

    pid_list.append(recv_hb_from_client_proc_hdl.pid)
    logger.info("Start sending my heartbeat to the client in every %d sec to %s"
                % (hearbeat_interval, WORKER_HB_Q_POSTFIX))

    # Start poison receive thread
    recv_poison_proc_hdl = mp.Process(target=recv_reproduce_or_die_proc,
                                      args=(pool_name_param,
                                            cluster_name,
                                            mem_per_node_to_request,
                                            mem_per_cpu_to_request,
                                            num_nodes_to_request,
                                            num_cpus_to_request,
                                            job_time_to_request,
                                            worker_clone_time_rate,
                                            num_workers_per_node,
                                            JTM_WORKER_POISON_EXCH,
                                            JTM_WORKER_POISON_Q))
    recv_poison_proc_hdl.daemon = True
    recv_poison_proc_hdl.start()
    pid_list.append(recv_poison_proc_hdl.pid)

    # Start hb send thread
    try:
        send_hb_to_client_proc_hdl = mp.Process(target=recv_hb_from_client_proc2,
                                                args=(inner_task_request_queue,
                                                      JTM_CLIENT_HB_EXCH,
                                                      CLIENT_HB_Q_POSTFIX))
        send_hb_to_client_proc_hdl.daemon = True
        send_hb_to_client_proc_hdl.start()
    except OSError as e:
        logger.exception("Worker termination request: {}".format(e))
        proc_clean()
        ch.stop_consuming()
        ch.close()
        conn.close()
        sys.exit(1)

    pid_list.append(send_hb_to_client_proc_hdl.pid)

    # Checking the total number of child processes
    try:
        check_processes_hdl = mp.Process(target=check_processes,
                                         args=(pid_list,))
        check_processes_hdl.daemon = True
        check_processes_hdl.start()
    except OSError as e:
        logger.exception("check_processes: {}".format(e))
        proc_clean()
        ch.stop_consuming()
        ch.close()
        conn.close()
        sys.exit(1)

    if worker_timeout_in_sec != 0:
        logger.info("The worker timeout is set to %s sec. Will not be terminated even without jtm's heartbeat.",
                    worker_timeout_in_sec)

    # Waiting for request
    ch.basic_qos(prefetch_count=1)
    try:
        ch.basic_consume(queue=inner_task_request_queue,
                         on_message_callback=on_request,
                         auto_ack=False)
    except OSError as err:
        logger.exception("Worker terminated: {}".format(err))
        proc_clean()
        ch.stop_consuming()
        ch.close()
        conn.close()
        sys.exit(1)

    # NEW
    # Ref) https://github.com/pika/pika/blob/1.0.1/examples/basic_consumer_threaded.py
    #      https://stackoverflow.com/questions/51752890/how-to-disable-heartbeats-with-pika-and-rabbitmq
    # threads = []
    # on_message_callback = functools.partial(on_task_request, args=(conn, threads))
    # ch.basic_consume(queue=inner_task_request_queue,
    #                  on_message_callback=on_message_callback)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        proc_clean()
        ch.stop_consuming()
        ch.close()
        conn.close()

    # Wait for all to complete
    # Note: prefetch_count=1 ==> #thread = 1
    # for thread in threads:
    #     thread.join()

    if ch:
        ch.close()
    if conn:
        conn.close()

    return 0
