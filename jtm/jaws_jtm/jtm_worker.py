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

    09.12.2018 3.1.0: Branched out "jtm"; Updated for pika=0.12.0; Changed "type" to "exchange_type"
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

    10.12.2018 1.4.0: Added recv_reproduce_or_die_thread; jtm_kill done
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
import multiprocessing
import shortuuid
import datetime
import argparse
import json
import threading
import subprocess
import socket
import time
import pika

from jaws_jtm.common import setup_custom_logger, logger
from jaws_jtm.config import PIKA_VER, \
    COMPUTE_RESOURCES, \
    VERSION, \
    JTM_WORKER_STREAM_LOGGING, \
    JTM_WORKER_FILE_LOGGING, \
    RMQ_PORT, \
    USER_NAME, \
    PRODUCTION, \
    JOBTIME, \
    CORI_CONSTRAINT, \
    CORI_KNL_CHARGE_ACCNT, \
    CORI_KNL_QOS, \
    CORI_CHARGE_ACCNT, \
    CORI_QOS, \
    WORKER_HB_SEND_INTERVAL, \
    WORKER_TIMEOUT, \
    JTM_INNER_REQUEST_Q, \
    CTR, \
    RMQ_HOST, \
    LAWRENCIUM_PARTITION, \
    JGI_CLOUD_PARTITION, \
    LAWRENCIUM_QOS, \
    JGI_CLOUD_QOS, \
    LAWRENCIUM_CHARGE_ACCNT, \
    JGI_CLOUD_ACCNT, \
    ENV_ACTIVATION, \
    JTM_INNER_MAIN_EXCH, \
    WORKER_HB_Q_POSTFIX, \
    JTM_LOG, \
    JTM_WORKER_POISON_Q, \
    JTM_WORKER_POISON_EXCH, \
    JTM_TASK_KILL_Q, \
    JTM_TASK_KILL_EXCH, \
    JOB_LOG, \
    TASK_TYPE, \
    CNAME, \
    CLIENT_HB_Q_POSTFIX, \
    JTM_CLIENT_HB_EXCH, \
    HB_MSG, \
    WORKER_TYPE, \
    JTM_HOST_NAME, \
    JTM_WORKER_HB_EXCH, \
    DONE_FLAGS, \
    FILE_CHECK_INT_INC, \
    FILE_CHECKING_MAX_TRIAL, \
    FILE_CHECK_INTERVAL
from jaws_jtm.lib.rabbitmqconnection import RmqConnectionHB
from jaws_jtm.lib.resourceusage import get_cpu_load, \
    get_runtime, get_pid_tree, get_virtual_memory_usage, \
    get_resident_memory_usage, get_total_mem_usage_per_node, \
    get_num_workers_on_node
from jaws_jtm.lib.run import make_dir, run_sh_command
from jaws_jtm.lib.msgcompress import zdumps, zloads


# This ipc pipe is to send task_id (by on_request())to send_hb_to_client_thread()
# when a task is requested
PIPE_TASK_ID_SEND, PIPE_TASK_ID_RECV = multiprocessing.Pipe()

print("JTM Worker, version: {}".format(VERSION))  # VERSION <- Config.py

# -------------------------------------------------------------------------------
# Globals
# -------------------------------------------------------------------------------
RECV_HB_FROM_CLIENT_PROC_HANDLE = None
SEND_HB_TO_CLIENT_PROC_HANDLE = None
RECV_POISON_PROC_HANDLE = None
TASK_KILL_PROC_HANDLE = None
USER_PROC_PROC_ID = multiprocessing.Value("i", 0)
ALLOW_EMPTY_OUTPUT_FILE = False  # False ==> if size(output file)==0, return error code
UNIQ_WORKER_ID = str(shortuuid.uuid())
IS_CLIENT_ALIVE = False
WORKER_START_TIME = datetime.datetime.now()
PARENT_PROCESS_ID = os.getpid()  # parent process id


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
    try:
        p = subprocess.Popen(user_task_cmd,
                             shell=True,
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
    stdout_str = p.communicate()[0]
    p.wait()
    p.poll()

    # Prepare result to send back
    return_msg["out_files"] = out_files
    return_msg["worker_id"] = UNIQ_WORKER_ID
    return_msg["host_name"] = socket.gethostname()

    if USER_PROC_PROC_ID.value == -9:
        return_msg["done_flag"] = DONE_FLAGS["failed with user termination"]
        return_msg["ret_msg"] = "Task cancelled."
    else:
        if p.returncode == 2:
            logger.info("input file not found.")
            return_msg["done_flag"] = "1"
            return_msg["ret_msg"] = "Input file not found."

        elif p.returncode == 0:
            logger.info("Task# %s completed!" % (str(task_id)))

            # Output file checking
            if len(out_files):
                ofs = out_files.split(",")
                logger.debug("Number of output files = %d.", len(ofs))
                out_file_list = []

                for i in range(len(ofs)):
                    out_file_list.append(ofs[i])

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
    # logger.debug("Result reply: %s" % str(json_data))
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

    # Send taskid to send_hb_to_client_thread() so that the task_id can be shown in the hb messages
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
        ch.stop_consuming()
        ch.close()

        if RECV_HB_FROM_CLIENT_PROC_HANDLE:
            RECV_HB_FROM_CLIENT_PROC_HANDLE.terminate()
        if SEND_HB_TO_CLIENT_PROC_HANDLE:
            SEND_HB_TO_CLIENT_PROC_HANDLE.terminate()
        if RECV_POISON_PROC_HANDLE:
            RECV_POISON_PROC_HANDLE.terminate()
        if TASK_KILL_PROC_HANDLE:
            TASK_KILL_PROC_HANDLE.terminate()

        sys.exit(0)

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
    thread = threading.Thread(target=run_something, args=(msg_unzipped, result_dict))
    thread.start()
    while thread.is_alive():  # Loop while the thread is processing
        ch._connection.sleep(1.0)
    json_data = json.dumps(result_dict)
    logger.debug("Result reply: %s" % str(json_data))
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

    # Send taskid=0 to send_hb_to_client_thread() b/c the requested task is completed
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

                # If "-z" option is used, allow zero sized output file
                if ALLOW_EMPTY_OUTPUT_FILE and file_size == 0:
                    file_size = 1

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
def send_hb_to_client_thread(root_proc_id, interval, queue_child_proc_id, worker_id, pipe_task_id,
                             slurm_job_id, worker_type, mem_per_node, mem_per_core, num_cores,
                             job_time, clone_time_rate, task_queue_name, pool_name, nwpn):
    """
    Send heartbeats to the client
    :param root_proc_id: root_proc_id of this worker
    :param interval: time interval to send heartbeats to the client
    :param queue_child_proc_id: child pid for user command
    :param worker_id: worker id
    :param pipe_task_id: Pythin IPC pipe for getting task id
    :param slurm_job_id: SLURM job id
    :param worker_type: worker type [manual, static, dynamic]
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

    exch_name = JTM_WORKER_HB_EXCH
    worker_hb_queue = WORKER_HB_Q_POSTFIX
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
            logger.debug("rootpid = {} childpid = {}".format(root_proc_id, queue_child_proc_id))
            root_proc_id = int(root_proc_id)
            if queue_child_proc_id.value == 0:
                child_pid = root_proc_id
            else:
                child_pid = int(queue_child_proc_id.value)

            # Collect pids from process tree
            root_pid_list = get_pid_tree(root_proc_id)
            logger.debug("get_pid_tree(root_proc_id={}) = {}".format(root_proc_id, root_pid_list))

            if root_proc_id != child_pid:
                if child_pid == -9:  # if it's terminated by "kill"
                    child_pid = root_proc_id
                pidListChild = get_pid_tree(child_pid)
                logger.debug("get_pid_tree(child_pid={}) = {}".format(child_pid, pidListChild))
                if len(pidListChild) > 0:
                    proc_id_list_merged = root_pid_list + pidListChild[1:]
                else:
                    # NOTE: be careful on this resetting! Might lose child pid
                    # queue_child_proc_id.value = 0
                    proc_id_list_merged = root_pid_list
            else:
                proc_id_list_merged = root_pid_list

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
            except Exception:
                logger.exception("get_cpu_load() exception")
                os._exit(1)

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
            num_workers_on_node = get_num_workers_on_node()

            # Check if there is any task id in the ipc pipe
            if pipe_task_id.poll():
                task_id = pipe_task_id.recv()

            today = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            end_date_time = None
            worker_age_in_minute = 0

            try:
                # Note: this only works with SLURM. For other HPCs/Clouds, methods to get the
                #  remaining life are needed! ==> resovled!
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
                    # worker_age_in_minute = -divmod(delta.total_seconds(), 60)[0]
                    # end_date_time = entries["EndTime"][0]

                    # NEW
                    # jobtime2 = datetime.datetime.strptime(job_time, '%H:%M:%S')
                    # end_date_time = WORKER_START_TIME + timedelta(seconds=jobtime2.second+jobtime2.minute*60+
                    # jobtime2.hour*3600)
                    # 24:00:00 --> seconds
                    job_runtime_in_sec = int(job_time.split(':')[0]) * 3600 + \
                                         int(job_time.split(':')[1]) * 60 + \
                                         int(job_time.split(':')[2])
                    end_date_time = WORKER_START_TIME + datetime.timedelta(seconds=job_runtime_in_sec)
                    delta = end_date_time - datetime.datetime.now()
                    worker_age_in_minute = divmod(delta.total_seconds(), 60)[0]
            except Exception as e:
                logger.critical("Something wrong in computing remaining wall clock time: %s", e)
                os._exit(1)

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

            ncore_param = multiprocessing.cpu_count()
            if WORKER_TYPE[worker_type] > 0:
                ncore_param = num_cores

            msg_dict_to_send = {HB_MSG["child_pid"]: child_pid,
                                HB_MSG["clone_time_rate"]: clone_time_rate,
                                HB_MSG["cpu_load"]: max_cpu_load,
                                HB_MSG["end_date"]: today,  # Note: for dynamic worker endDate update
                                HB_MSG["host_name"]: host_name,
                                HB_MSG["ip_address"]: ip_address,
                                HB_MSG["job_time"]: job_time if WORKER_TYPE[worker_type] > 0 else None,
                                HB_MSG["jtm_host_name"]: JTM_HOST_NAME,
                                HB_MSG["life_left"]: worker_age_in_minute,
                                HB_MSG["mem_per_core"]: mem_per_core if WORKER_TYPE[worker_type] > 0 else "",
                                HB_MSG["mem_per_node"]: mem_per_node if WORKER_TYPE[worker_type] > 0 else "",
                                HB_MSG["num_cores"]: ncore_param,
                                HB_MSG["num_tasks"]: hb_queue_len,
                                HB_MSG["num_workers_on_node"]: num_workers_on_node,
                                HB_MSG["perc_mem_used"]: perc_used_mem,
                                HB_MSG["pool_name"]: pool_name,
                                HB_MSG["ret_msg"]: "hb",
                                HB_MSG["rmem_usage"]: rmem_usage,
                                HB_MSG["root_pid"]: root_proc_id,
                                HB_MSG["run_time"]: proc_run_time,
                                HB_MSG["slurm_jobid"]: slurm_job_id,
                                HB_MSG["task_id"]: task_id,
                                HB_MSG["vmem_usage"]: vmem_usage,
                                HB_MSG["worker_id"]: worker_id,
                                HB_MSG["worker_type"]: WORKER_TYPE[worker_type],
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
            os._exit(1)

        # Todo: Need to start with shorted interval like (0.5-1 sec) for about 30 sec
        #  from the beginning and then use the user specified interval
        # time.sleep(float(interval))
        conn.process_data_events(time_limit=float(interval))

    # unreachable
    ch.close()
    conn.close()


# -----------------------------------------------------------------------
def recv_hb_from_client_thread2(task_queue, worker_id):
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()

    exch_name = JTM_CLIENT_HB_EXCH
    # hb_interval = WORKER_HB_RECV_INTERVAL
    host_name = socket.gethostname()

    # Declare exchange
    ch.exchange_declare(exchange=exch_name,
                        # exchange_type="fanout",
                        exchange_type="topic",
                        durable=False,
                        auto_delete=False)

    # Declare queue
    client_hb_queue_name = "_jtm_worker_%s%s" % (worker_id, CLIENT_HB_Q_POSTFIX)
    ch.queue_declare(queue=client_hb_queue_name,
                     durable=False,
                     exclusive=True,
                     auto_delete=True)

    # Bind exchange to hb queue for receiving the client hb
    ch.queue_bind(exchange=exch_name,
                  queue=client_hb_queue_name,
                  routing_key="*." + CNAME)

    logger.info('[*] Waiting for the manager.')

    def stop():
        logger.info("Terminate myself: {} on {}".format(worker_id, host_name))
        if RECV_HB_FROM_CLIENT_PROC_HANDLE:
            RECV_HB_FROM_CLIENT_PROC_HANDLE.terminate()
        if SEND_HB_TO_CLIENT_PROC_HANDLE:
            SEND_HB_TO_CLIENT_PROC_HANDLE.terminate()
        if RECV_POISON_PROC_HANDLE:
            RECV_POISON_PROC_HANDLE.terminate()
        if TASK_KILL_PROC_HANDLE:
            TASK_KILL_PROC_HANDLE.terminate()
        ch.stop_consuming()
        # os.kill(int(ppid), signal.SIGTERM)  # term jtm-worker process
        # sys.exit(10)
        os._exit(1)

    def callback(ch, method, properties, body):
        msg_unzipped = json.loads(zloads(body))
        if msg_unzipped["task_type"] == TASK_TYPE["term"] and msg_unzipped["task_queue"] == task_queue:
            stop()
        elif msg_unzipped["task_type"] == TASK_TYPE["hb"]:
            # Poison from the manager. Terminate itself.
            logger.debug("Received heartbeats from the client, %r:%r" % (method.routing_key, msg_unzipped))

    ch.basic_consume(queue=client_hb_queue_name,
                     on_message_callback=callback,
                     auto_ack=True)

    ch.start_consuming()


# -------------------------------------------------------------------------------
def recv_task_kill_request_thread(task_queue, worker_id, cluster_name):
    """
    Wait for task termination request from JTM
    :param task_queue: task queue (pool) name
    :param worker_id: worker id
    :param cluster_name: cluster name
    :return:
    """
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()

    exch_name = JTM_TASK_KILL_EXCH
    queue_name = JTM_TASK_KILL_Q
    routing_key = str(worker_id)

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

        if msg_unzipped["worker_id"] == worker_id:
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
                    so, se, ec = run_sh_command(kill_cmd, live=True, log=logger)
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
                so, se, ec = run_sh_command(kill_cmd, live=True, log=logger)
                if ec == 0:
                    logger.info("Successfully terminate a user task process.")
                else:
                    logger.warning("User process not found. Failed to execute the command, %s" % (kill_cmd))

            else:
                logger.warning("No valid child process id to terminate.")

            ch.basic_ack(delivery_tag=method.delivery_tag)

        else:
            # ch.basic_nack(delivery_tag=method.delivery_tag)
            # If worker_id is not for me, reject and requeue it
            ch.basic_reject(delivery_tag=method.delivery_tag, requeue=True)

    # Waiting for a task kill
    ch.basic_qos(prefetch_count=1)

    if int(PIKA_VER[0]) < 1:  # v0.13.1
        ch.basic_consume(on_kill, queue=queue_name, no_ack=False)
    else:  # v1.0.1 or higher
        ch.basic_consume(queue=queue_name, on_message_callback=on_kill, auto_ack=False)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        ch.stop_consuming()


# -------------------------------------------------------------------------------
def recv_reproduce_or_die_thread(task_queue, worker_id, cluster_name, mem_per_node, mem_per_core,
                                 num_nodes, num_cores, wallclocktime, clone_time_rate, nwpn):
    """
    Wait for process termination request from JTM
    :param task_queue: task queue (pool) name
    :param worker_id: worker id
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

    exch_name = JTM_WORKER_POISON_EXCH
    queue_name = JTM_WORKER_POISON_Q

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
                  routing_key=worker_id)

    def on_poison(ch, method, props, body):
        msg_unzipped = json.loads(zloads(body))

        if msg_unzipped["worker_id"] == worker_id:

            # If static worker cloning request
            if msg_unzipped["num_clones"] == -1:
                # Kill the worker request
                pass
            elif msg_unzipped["num_clones"] == 1:
                logger.info("Static worker cloning request received.")

                jtm_worker_cmd = "{} && jtm-worker -wt static -cl {} -t {} -ct {}".format(
                                ENV_ACTIVATION, cluster_name, wallclocktime, clone_time_rate)

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
                _, _, _ = run_sh_command(jtm_worker_cmd + " --dry-run", live=True, log=logger)

            # Note: 01072019 no more auto cloning of dynamic worker
            # elif msg_unzipped["num_clones"] > 1:
            #     logger.info("Dynamic worker cloning request received.")
            #     jtm_worker_cmd = "{} && jtm-worker -wt dynamic -cl {} -t {} -ct {}".format(
            #                     ENV_ACTIVATION, cluster_name, wallclocktime, clone_time_rate)
            #     if task_queue is not None:
            #         jtm_worker_cmd += " -p {}".format(task_queue)
            #     if num_nodes != "" and int(num_nodes) > 0:
            #         jtm_worker_cmd += " -N {} -m {} -c {} -nw {}".format(num_nodes, mem_per_node, num_cores, nwpn)
            #     else:
            #         jtm_worker_cmd += " -c {}".format(num_cores)
            #         if mem_per_node not in ("", None):  # if node mem defined, set --mem
            #             jtm_worker_cmd += " -m {}".format(mem_per_node)
            #         else:
            #             jtm_worker_cmd += " -mc {}".format(mem_per_core)
            #
            #     for _ in range(0, msg_unzipped["num_clones"]):
            #         logger.info("Executing {}".format(jtm_worker_cmd))
            #         _, _, _ = run_sh_command(jtm_worker_cmd, live=True, log=logger)
            #         # _, _, _ = run_sh_command(jtm_worker_cmd + " --dry-run", live=True, log=logger)

            ch.basic_ack(delivery_tag=method.delivery_tag)

        else:
            # ch.basic_nack(delivery_tag=method.delivery_tag)
            # If worker_id is not for me, reject and requeue it
            ch.basic_reject(delivery_tag=method.delivery_tag, requeue=True)

    # waiting for poison or cloning command
    ch.basic_qos(prefetch_count=1)
    if int(PIKA_VER[0]) < 1:  # v0.13.1
        ch.basic_consume(on_poison, queue=queue_name, no_ack=False)
    else:  # v1.0.1 or higher
        ch.basic_consume(queue=queue_name, on_message_callback=on_poison, auto_ack=False)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        ch.stop_consuming()


# -------------------------------------------------------------------------------
def worker():
    desc = u"JTM worker"
    parser = argparse.ArgumentParser(description=desc)
    # parser.add_argument("-cj", "--cromwell-jid",
    #                     help="Set Cromwell job name",
    #                     dest="cromwell_job_id")
    parser.add_argument("-hb", "--heartbeat",
                        help="heartbeat interval in second (default=10).",
                        dest="hearbeat_interval",
                        type=int)  # default 5 sec
    parser.add_argument("-l", "--loglevel",
                        help="Set loglevel (default=info).",
                        dest="log_level",
                        default="info")
    parser.add_argument("-jd", "--job_dir",
                        help="jtm job file path.",
                        dest="job_script_dir_name",
                        default=JOB_LOG)
    parser.add_argument("-ld", "--logdir",
                        help="jtm log file path.",
                        dest="log_dir_name",
                        default=JTM_LOG)
    parser.add_argument("-p", "--pool-name",
                        help="Set user-defined pool name. This should be same with task's 'pool' name.",
                        dest="task_queue_name",
                        default="small")
    parser.add_argument("-to", "--timeout",
                        help="Set the timer for worker to terminate. If there is no request from the client \
                        for the specified seconds, the worker terminates itself (default: 60 seconds)",
                        dest="worker_timeout_in_sec",
                        type=int)
    parser.add_argument("-zf", "--zerofile",
                        action="store_true",
                        help="Allow zero size output file(s) (default=False).",
                        dest="zero_size_file",
                        default=False)
    parser.add_argument("-dr", "--dry-run",
                        action="store_true",
                        help="Dry run. To print batch job script on screen",
                        dest="dry_run",
                        default=False)
    parser.add_argument("-j", "--jobid",
                        help="Slurm job id",
                        dest="slurm_job_id",
                        type=int,
                        default=0)

    # worker type params
    # -------------------
    parser.add_argument("-wt", "--worker-type",
                        help="jtm worker type = [static | dynamic]",
                        dest="worker_type",
                        default="manual",
                        choices=["static", "dynamic"])
    parser.add_argument("-cl", "--cluster",
                        help="cluster name = [%s]" % (str(COMPUTE_RESOURCES)),
                        dest="cluster_name",
                        default="local",
                        choices=COMPUTE_RESOURCES)
    parser.add_argument("-ct", "--clone-time",
                        help="cloning timing",
                        dest="worker_clone_time_rate",
                        type=float)
    parser.add_argument("-nw", "--num-workers",
                        help="number of worker per node. Default=1",
                        dest="num_workers_per_node",
                        type=int)
    parser.add_argument("-wi", "--worker-id",
                        help="worker id for dynamic workers",
                        dest="worker_id")

    # slurm related params
    # -------------------
    parser.add_argument("-A", "--account",
                        help="Charging account.",
                        dest="charging_account")
    parser.add_argument("-N", "--nodes",
                        help="number of nodes",
                        dest="num_nodes_to_request",
                        type=int)
    parser.add_argument("-c", "--cpus-per-task",
                        help="number of cpus",
                        dest="num_cores_to_request",
                        type=int)
    parser.add_argument("-C", "--constraint",
                        help="Cori constraint, haswell or knl (default: haswell)",
                        dest="constraint",
                        choices=["haswell", "knl", "skylake"],
                        default="haswell")
    parser.add_argument("-m", "--mem",
                        help="specify the real memory required per node",
                        dest="mem_per_node_to_request")
    parser.add_argument("-mc", "--mem-per-cpu",
                        help="minimum memory required per allocated CPU",
                        dest="mem_per_cpu_to_request")
    parser.add_argument("-q", "--qos",
                        help="quality of service",
                        dest="qos")
    parser.add_argument("-t", "--time",
                        help="limit on the total run time",
                        dest="job_time_to_request")
    parser.add_argument("-J", "--job-name",
                        help="Slurm job name",
                        dest="job_name")

    parser.add_argument("-v", "--version",
                        action="version",
                        version=VERSION)  # VERSION <- Config.py
    args = parser.parse_args()

    # Note: this is for allowing zero sized output file when check if the output file is
    #  successfully created. Currently output file checking is opted out in JTM.
    global ALLOW_EMPTY_OUTPUT_FILE
    ALLOW_EMPTY_OUTPUT_FILE = args.zero_size_file

    # Set uniq worker id if worker id is provided in the params
    if args.worker_id:
        global UNIQ_WORKER_ID
        UNIQ_WORKER_ID = args.worker_id

    # Logger setting
    log_dir_name = args.log_dir_name
    # make_dir_p(log_dir_name)
    make_dir(log_dir_name)

    # Job dir setting
    job_script_dir_name = args.job_script_dir_name
    # make_dir_p(job_script_dir_name)
    make_dir(job_script_dir_name)

    setup_custom_logger(args.log_level, log_dir_name,
                        JTM_WORKER_STREAM_LOGGING, JTM_WORKER_FILE_LOGGING,
                        worker_id=UNIQ_WORKER_ID)

    logger.info("Set jtm log file location to %s", args.log_dir_name)
    logger.info("Set jtm job file location to %s", args.job_script_dir_name)
    logger.info("RabbitMQ broker: %s", RMQ_HOST)
    logger.info("RabbitMQ port: %s", RMQ_PORT)
    logger.info("Pika version: %s", PIKA_VER)
    logger.info("JTM user name: %s", USER_NAME)
    logger.info("Unique worker ID: %s", UNIQ_WORKER_ID)
    logger.info("**************")
    logger.info("Run mode: %s", "prod" if PRODUCTION else "dev")
    logger.info("**************")

    # Slurm config
    num_nodes_to_request = 0
    if args.num_nodes_to_request:
        num_nodes_to_request = args.num_nodes_to_request
        # Todo
        #  Cori and JGI Cloud are exclusive allocation. So this is not needed.
        # assert args.mem_per_node_to_request is not None, "-N needs --mem-per-cpu (-mc) setting."

    ###########################################################################
    # 11.13.2018 decided to remove all default values from argparse
    num_workers_per_node = args.num_workers_per_node if args.num_workers_per_node else 1
    num_cpus_to_reques = args.num_cores_to_request if args.num_cores_to_request else 1
    mem_per_cpu_to_request = args.mem_per_cpu_to_request if args.mem_per_cpu_to_request else "1GB"
    mem_per_node_to_request = args.mem_per_node_to_request if args.mem_per_node_to_request else ""
    ###########################################################################

    job_time_to_request = args.job_time_to_request if args.job_time_to_request else JOBTIME

    constraint = args.constraint if args.constraint else CORI_CONSTRAINT
    if constraint == "knl":
        charging_account = args.charging_account if args.charging_account else CORI_KNL_CHARGE_ACCNT
        qos = args.qos if args.qos else CORI_KNL_QOS
    else:
        charging_account = args.charging_account if args.charging_account else CORI_CHARGE_ACCNT
        qos = args.qos if args.qos else CORI_QOS

    worker_type = args.worker_type
    job_name = "jtm_worker_" + args.task_queue_name

    # Set task queue name
    inner_task_request_queue = None
    hearbeat_interval = WORKER_HB_SEND_INTERVAL
    if args.hearbeat_interval:
        hearbeat_interval = args.hearbeat_interval
    worker_timeout_in_sec = WORKER_TIMEOUT
    if args.worker_timeout_in_sec:
        worker_timeout_in_sec = args.worker_timeout_in_sec

    # If you want to create a custom worker pool, use this.
    # JTM detects "pool" field from task json, and creates custom pool if it's set.
    # All tasks with the pool name will be directed to the custom pool
    # user_account_name = getpass.getuser().lower()

    # Todo: feature to add a manual worker to a specific cluster??
    #  ex) $ jtm-worker -cl cori -p small (on cori)
    # if args.worker_type == "manual" and args.cluster_name != "local":
    #     # inner_task_request_queue = "_jtm_inner_request_queue." + args.cluster_name + "." + args.task_queue_name
    #     JTM_INNER_REQUEST_Q = "_jtm_inner_request_queue." + args.cluster_name + "." + USER_NAME

    if args.task_queue_name:
        inner_task_request_queue = JTM_INNER_REQUEST_Q + "." + args.task_queue_name
    else:  # not reachable. default inner_task_request_queue = small
        inner_task_request_queue = JTM_INNER_REQUEST_Q

    worker_clone_time_rate = args.worker_clone_time_rate if args.worker_clone_time_rate else CTR
    if worker_type in ("static", "dynamic"):
        assert args.cluster_name != "" and \
               args.cluster_name != "local", "Static or dynamic worker needs a cluster setting (-cl)."

    slurm_job_id = args.slurm_job_id
    cluster_name = args.cluster_name  # default=local
    dry_run = args.dry_run

    if cluster_name == "cori" and mem_per_cpu_to_request != "" and \
            float(mem_per_cpu_to_request.replace("GB", "").replace("G", "")) > 1.0:
        logger.critical("--mem-per-cpu in Cori shouldn't be larger than 1GB. User '--mem' instead.")
        sys.exit(1)

    logger.info("RabbitMQ broker: %s", RMQ_HOST)
    logger.info("Task queue name: %s", inner_task_request_queue)
    logger.info("Worker type: %s", args.worker_type)

    # sbatch static worker
    # -N: If -N is not specified, the default behavior is to allocate enough nodes to satisfy
    #     the requirements of the -n and -c options. Recommended to set --mem
    # -c: Set the num cores needed. Recommended to set --mem-per-cpu

    # SBATCH
    # -N, --nodes=<minnodes[-maxnodes]>
    #               Request  that  a minimum of minnodes nodes be allocated to this job.  A maximum node count may also
    #               be specified with maxnodes.  If only one number is specified, this is used as both the minimum  and
    #               maximum  node count.  The partition's node limits supersede those of the job.  If a job's node lim-
    #               its are outside of the range permitted for its associated partition, the job  will  be  left  in  a
    #               PENDING  state.   This  permits  possible  execution  at  a later time, when the partition limit is
    #               changed.  If a job node limit exceeds the number of nodes configured in the partition, the job will
    #               be  rejected.  Note that the environment variable SLURM_JOB_NODES will be set to the count of nodes
    #               actually allocated to the job. See the ENVIRONMENT VARIABLES  section for more information.  If  -N
    #               is  not  specified, the default behavior is to allocate enough nodes to satisfy the requirements of
    #               the -n and -c options.  The job will be allocated as many nodes as possible within the range speci-
    #               fied  and  without  delaying the initiation of the job.  The node count specification may include a
    #               numeric value followed by a suffix of "k" (multiplies numeric value by 1,024)  or  "m"  (multiplies
    #               numeric value by 1,048,576).
    # -c, --cpus-per-task=<ncpus>
    #               Advise  the  Slurm  controller  that  ensuing job steps will require ncpus number of processors per
    #               task.  Without this option, the controller will just try to allocate one processor per task.
    #
    #               For instance, consider an application that has 4 tasks, each requiring 3 processors.  If our  clus-
    #               ter is comprised of quad-processors nodes and we simply ask for 12 processors, the controller might
    #               give us only 3 nodes.  However, by using the --cpus-per-task=3 options, the controller  knows  that
    #               each  task requires 3 processors on the same node, and the controller will grant an allocation of 4
    #               nodes, one for each of the 4 tasks.
    # -n, --ntasks=<number>
    #               sbatch  does  not  launch tasks, it requests an allocation of resources and submits a batch script.
    #               This option advises the Slurm controller that job steps run within the  allocation  will  launch  a
    #               maximum of number tasks and to provide for sufficient resources.  THE DEFAULT IS ONE TASK PER NODE,
    #               BUT NOTE THAT THE --cpus-per-task OPTION WILL CHANGE THIS DEFAULT.
    # --mem=<size[units]>
    #               Specify the real memory required per node.  Default units are megabytes unless the SchedulerParame-
    #               ters  configuration  parameter includes the "default_gbytes" option for gigabytes.  Different units
    #               can be specified using the suffix [K|M|G|T].  Default value is DefMemPerNode and the maximum  value
    #               is  MaxMemPerNode.  If  configured, both parameters can be seen using the scontrol show config com-
    #               mand.  THIS PARAMETER WOULD GENERALLY BE USED  IF  WHOLE  NODES  ARE  ALLOCATED  TO  JOBS  (Select-
    #               Type=select/linear).  Also see --mem-per-cpu.  --mem and --mem-per-cpu are mutually exclusive.
    #
    #               NOTE: A memory size specification of zero is treated as a special case and grants the job access to
    #               all of the memory on each node.  If the job is allocated multiple nodes in a heterogeneous cluster,
    #               the  memory  limit on each node will be that of the node in the allocation with the smallest memory
    #               size (same limit will apply to every node in the job's allocation).
    #
    #               NOTE: Enforcement of memory limits currently relies upon the task/cgroup plugin or enabling of
    #               accounting, which samples memory use on a periodic basis (data need not be stored, just collected).
    #               In both cases memory use is based upon the job's Resident Set Size (RSS). A  task  may  exceed  the
    #               memory limit until the next periodic accounting sample.
    # --mem-per-cpu=<size[units]>
    #               Minimum memory required per allocated CPU.  Default units are megabytes unless the SchedulerParame-
    #               ters configuration parameter includes the "default_gbytes" option for gigabytes.  Default value  is
    #               DefMemPerCPU  and  the  maximum  value  is  MaxMemPerCPU (see exception below). If configured, both
    #               parameters can  be  seen  using  the  scontrol  show  config  command.  Note  that  if  the  job's
    #               --mem-per-cpu value exceeds the configured MaxMemPerCPU, then the user's limit will be treated as a
    #               memory limit per task; --mem-per-cpu will be reduced  to  a  value  no  larger  than  MaxMemPerCPU;
    #               --cpus-per-task  will  be  set and the value of --cpus-per-task multiplied by the new --mem-per-cpu
    #               value will equal the original --mem-per-cpu value specified by the user.  This parameter would gen-
    #               erally  be  used  if  individual processors are allocated to jobs (SelectType=select/cons_res).  If
    #               resources are allocated by the core, socket or whole nodes; the number of CPUs allocated to  a  job
    #               may  be  higher  than the task count and the value of --mem-per-cpu should be adjusted accordingly.
    #               Also see --mem.  --mem AND --mem-per-cpu ARE MUTUALLY EXCLUSIVE.
    #
    # NOTE: Cori shared queue doesn't need "account" and "qos"
    #       If you set --qos=genepool or --qos=genepool_shared, you don't need to set "-C haswell"
    #                                                           you should set "-A fungalp"
    #       If you set -A m342, no qos needed but you need to set "-C haswell"
    #
    if slurm_job_id == 0 and worker_type in ["static", "dynamic"]:
        batch_job_script_file = os.path.join(job_script_dir_name, "jtm_%s_worker_%s.job" %
                                             (worker_type, UNIQ_WORKER_ID))
        batch_job_script_str = ""
        batch_job_misc_params = ""

        if cluster_name in ("cori", "lawrencium", "jgi_cloud", "jaws_lbl_gov"):

            with open(batch_job_script_file, "w") as jf:
                batch_job_script_str += "#!/bin/bash -l"

                if cluster_name in ("cori"):

                    if args.num_nodes_to_request:
                        batch_job_script_str += """
#SBATCH -N %(num_nodes_to_request)d
#SBATCH --mem=%(mem)s""" % dict(num_nodes_to_request=num_nodes_to_request, mem=mem_per_node_to_request)
                        batch_job_misc_params += " -N %(num_nodes_to_request)d -m %(mem)s" % \
                                                 dict(num_nodes_to_request=num_nodes_to_request,
                                                      mem=mem_per_node_to_request)

                        if args.num_cores_to_request:
                            batch_job_script_str += """
#SBATCH -c %(num_cores)d""" % dict(num_cores=num_cpus_to_reques)
                            batch_job_misc_params += " -c %(num_cores)d" % \
                                                     dict(num_cores=num_cpus_to_reques)

                    else:
                        batch_job_script_str += """
#SBATCH -c %(num_cores)d""" % dict(num_cores=num_cpus_to_reques)
                        batch_job_misc_params += " -c %(num_cores)d" % \
                                                 dict(num_cores=num_cpus_to_reques)

                        if args.mem_per_node_to_request:
                            batch_job_script_str += """
#SBATCH --mem=%(mem)s""" % dict(mem=mem_per_node_to_request)
                            batch_job_misc_params += " -m %(mem)s " % \
                                                     dict(mem=mem_per_node_to_request)
                        else:
                            batch_job_script_str += """
#SBATCH --mem-per-cpu=%(mempercore)s""" % dict(mempercore=mem_per_cpu_to_request)
                            batch_job_misc_params += " -mc %(mempercore)s" % \
                                                     dict(mempercore=mem_per_cpu_to_request)

                        if args.worker_id:
                            batch_job_misc_params += " -wi %(worker_id)s${i}" % \
                                                     dict(worker_id=UNIQ_WORKER_ID)

                    ###########################
                    if cluster_name == "cori":

                        # Need to set both --qos=genepool (or genepool_shared) _and_ -A fungalp
                        # OR
                        # no qos _and_ -A m342 _and_ -C haswell

                        # Note: currently constraint in ["haswell" | "knl"]
                        if args.constraint == "haswell":
                            if args.qos:
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

                        elif args.constraint == "knl":
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

                        elif args.constraint == "skylake":
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
                        if args.constraint != "skylake":
                            excl_param = "#SBATCH --exclusive"

                        tq_param = ""
                        if args.task_queue_name:
                            tq_param = "-p " + args.task_queue_name

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
    jtm-worker --jobid $SLURM_JOB_ID \
-cl cori \
-wt %(worker_type)s \
-t %(wall_time)s \
-ct %(clone_time_rate)f %(task_queue)s \
-nw %(num_workers_per_node)d \
-C %(constraint)s \
%(other_params)s &
    sleep 1
done
wait
""" % \
                                                dict(wall_time=job_time_to_request,
                                                     job_dir=job_script_dir_name,
                                                     worker_id=UNIQ_WORKER_ID,
                                                     worker_type=worker_type,
                                                     clone_time_rate=worker_clone_time_rate,
                                                     task_queue=tq_param,
                                                     num_workers_per_node=num_workers_per_node,
                                                     env_activation_cmd=ENV_ACTIVATION,
                                                     other_params=batch_job_misc_params,
                                                     constraint=args.constraint,
                                                     job_name=job_name,
                                                     exclusive=excl_param)

                elif cluster_name in ("lawrencium", "jgi_cloud", "jaws_lbl_gov"):

                    if args.worker_id:
                        # batch_job_misc_params += " -wi %s" % (args.worker_id)
                        batch_job_misc_params += " -wi %(worker_id)s${i}" % dict(worker_id=UNIQ_WORKER_ID)

                    tp_param = ""
                    if args.task_queue_name:
                        tp_param = "-p " + args.task_queue_name
                    part_param = ""
                    if cluster_name == "lawrencium":
                        part_param = LAWRENCIUM_PARTITION
                    else:
                        part_param = JGI_CLOUD_PARTITION
                    qos_param = ""
                    if cluster_name == "lawrencium":
                        qos_param = LAWRENCIUM_QOS
                    else:
                        qos_param = JGI_CLOUD_QOS
                    charge_param = ""
                    if cluster_name == "lawrencium":
                        charge_param = LAWRENCIUM_CHARGE_ACCNT
                    else:
                        charge_param = JGI_CLOUD_ACCNT
                    nnode_param = 1
                    if args.num_nodes_to_request:
                        nnode_param = num_nodes_to_request

                    mnode_param = "#SBATCH --mem=%(mem)s" % dict(mem=mem_per_node_to_request)

                    setenv_param = ""
                    if cluster_name == "jaws_lbl_gov":
                        setenv_param = "export JTM_HOST_NAME=jaws_lbl_gov"

                    batch_job_script_str += """
#SBATCH --time=%(wall_time)s
#SBATCH --job-name=%(job_name)s
#SBATCH --partition=%(partitio_name)s
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
    jtm-worker --jobid $SLURM_JOB_ID \
-cl %(lbl_cluster_name)s \
-wt %(worker_type)s \
-t %(wall_time)s \
-ct %(clone_time_rate)f %(task_queue)s \
-nw %(num_workers_per_node)d \
%(other_params)s &
    sleep 1
done
wait
""" % \
                                            dict(wall_time=job_time_to_request,
                                                 job_name=job_name,
                                                 partitio_name=part_param,
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
                                                 worker_type=worker_type,
                                                 clone_time_rate=worker_clone_time_rate,
                                                 task_queue=tp_param,
                                                 other_params=batch_job_misc_params)

                jf.writelines(batch_job_script_str)

            os.chmod(batch_job_script_file, 0o775)

            if dry_run:
                print(batch_job_script_str)
                sys.exit(0)

            so, se, ec = run_sh_command("sbatch --parsable %s" % (batch_job_script_file), live=True, log=logger)
            assert ec == 0, "Failed to run jtm-worker to sbatch dynamic worker."
            return ec
            # if ec == 0:
            #     sys.exit(0)
            # else:
            #     logger.critical("Failed to run jtm-worker to sbatch dynamic worker.")
            #     sys.exit(1)

        elif cluster_name == "aws":
            pass

    # If it's spawned by sbatch
    # Todo: need to record job_id, worker_id, worker_type, starting_time, wallclocktime
    # scontrol show jobid -dd <jobid> ==> EndTime
    # scontrol show jobid <jobid> ==> EndTime
    # sstat --format=AveCPU,AvePages,AveRSS,AveVMSize,JobID -j <jobid> --allsteps
    #
    # if endtime - starttime <= 10%, execute sbatch again
    # if slurm_job_id != 0 and worker_type == "static":
    #     logger.debug("worker_type: {}".format(worker_type))
    #     logger.debug("slurm_job_id: {}".format(slurm_job_id))

    # Dynamic workers creates [[two]] children when it approaches to the wallclocktime limit
    # considering the task queue length
    # Also, maintain the already requested number of workers
    # if no more workers needed, it won't call sbatch
    # elif slurm_job_id != 0 and worker_type == "dynamic":
    #     logger.debug("worker_type: {}".format(worker_type))
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

    # Start task termination thread
    global TASK_KILL_PROC_HANDLE
    TASK_KILL_PROC_HANDLE = multiprocessing.Process(target=recv_task_kill_request_thread,
                                                    args=(args.task_queue_name,
                                                          UNIQ_WORKER_ID,
                                                          cluster_name))
    TASK_KILL_PROC_HANDLE.daemon = True
    TASK_KILL_PROC_HANDLE.start()

    # Start hb receive thread
    global RECV_HB_FROM_CLIENT_PROC_HANDLE
    tp_name = ""
    if args.task_queue_name:
        tp_name = args.task_queue_name
    RECV_HB_FROM_CLIENT_PROC_HANDLE = multiprocessing.Process(target=send_hb_to_client_thread,
                                                              args=(PARENT_PROCESS_ID,
                                                                    hearbeat_interval,
                                                                    USER_PROC_PROC_ID,
                                                                    UNIQ_WORKER_ID,
                                                                    PIPE_TASK_ID_RECV,
                                                                    slurm_job_id,
                                                                    worker_type,
                                                                    mem_per_node_to_request,
                                                                    mem_per_cpu_to_request,
                                                                    num_cpus_to_reques,
                                                                    job_time_to_request,
                                                                    worker_clone_time_rate,
                                                                    inner_task_request_queue,
                                                                    tp_name,
                                                                    num_workers_per_node))

    RECV_HB_FROM_CLIENT_PROC_HANDLE.daemon = True
    RECV_HB_FROM_CLIENT_PROC_HANDLE.start()
    logger.info("Start sending my heartbeat to the client in every %d sec to %s" %
                (hearbeat_interval, WORKER_HB_Q_POSTFIX))

    # Start poison receive thread
    global RECV_POISON_PROC_HANDLE
    RECV_POISON_PROC_HANDLE = multiprocessing.Process(target=recv_reproduce_or_die_thread,
                                                      args=(args.task_queue_name,
                                                            UNIQ_WORKER_ID,
                                                            cluster_name,
                                                            mem_per_node_to_request,
                                                            mem_per_cpu_to_request,
                                                            num_nodes_to_request,
                                                            num_cpus_to_reques,
                                                            job_time_to_request,
                                                            worker_clone_time_rate,
                                                            num_workers_per_node))
    RECV_POISON_PROC_HANDLE.daemon = True
    RECV_POISON_PROC_HANDLE.start()

    # Start hb send thread
    global SEND_HB_TO_CLIENT_PROC_HANDLE
    SEND_HB_TO_CLIENT_PROC_HANDLE = multiprocessing.Process(target=recv_hb_from_client_thread2,
                                                            args=(inner_task_request_queue,
                                                                  UNIQ_WORKER_ID))
    SEND_HB_TO_CLIENT_PROC_HANDLE.daemon = True
    SEND_HB_TO_CLIENT_PROC_HANDLE.start()

    if worker_timeout_in_sec != 0:
        logger.info("The worker timeout is set to %s sec. Will not be terminated even without jtm's heartbeat.",
                    worker_timeout_in_sec)

    # Waiting for request
    ch.basic_qos(prefetch_count=1)

    # OLD
    if int(PIKA_VER[0]) < 1:  # v0.13.1
        ch.basic_consume(on_request, queue=inner_task_request_queue, no_ack=False)
    else:  # v1.0.1 or higher
        ch.basic_consume(queue=inner_task_request_queue, on_message_callback=on_request, auto_ack=False)

    # NEW
    # Ref) https://github.com/pika/pika/blob/1.0.1/examples/basic_consumer_threaded.py
    #      https://stackoverflow.com/questions/51752890/how-to-disable-heartbeats-with-pika-and-rabbitmq
    # threads = []
    # on_message_callback = functools.partial(on_task_request, args=(conn, threads))
    # if int(PIKA_VER[0]) < 1:  # v0.13.1
    #     ch.basic_consume(inner_task_request_queue, on_message_callback)
    # else:  # v1.0.1 or higher
    #     ch.basic_consume(queue=inner_task_request_queue, on_message_callback=on_message_callback)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        if RECV_HB_FROM_CLIENT_PROC_HANDLE:
            RECV_HB_FROM_CLIENT_PROC_HANDLE.terminate()
        if SEND_HB_TO_CLIENT_PROC_HANDLE:
            SEND_HB_TO_CLIENT_PROC_HANDLE.terminate()
        if RECV_POISON_PROC_HANDLE:
            RECV_POISON_PROC_HANDLE.terminate()
        if TASK_KILL_PROC_HANDLE:
            TASK_KILL_PROC_HANDLE.terminate()
        ch.stop_consuming()

    # Wait for all to complete
    # Note: prefetch_count=1 ==> #thread = 1
    # for thread in threads:
    #     thread.join()

    # Unreachable
    if ch:
        ch.close()
    if conn:
        conn.close()

    return 0


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    sys.exit(worker(sys.argv))
