#! /usr/bin/env python
# pylint: disable=C0111,C0103,R0205
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

    06.11.2019 5.4.0: Run run_user_task() with threading.Thread; Added ch._connection.sleep() to fix
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
    :return:
    """
    # Uncompress msg to get a taska
    logger.info(msg_unzipped)
    task_id = msg_unzipped["task_id"]
    user_task_cmd = msg_unzipped["user_cmd"]
    out_files = msg_unzipped["output_files"]
    task_type = msg_unzipped["task_type"]

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

    p = None
    time_out_in_minute = 0

    if THIS_WORKER_TYPE != "manual":
        # wait until WORKER_LIFE_LEFT_IN_MINUTE is updated
        if WORKER_LIFE_LEFT_IN_MINUTE.value <= 0:
            while True:
                logger.debug("worker life: %d", WORKER_LIFE_LEFT_IN_MINUTE.value)
                if WORKER_LIFE_LEFT_IN_MINUTE.value > 0:
                    break
                # time.sleep(1)
                ch._connection.process_data_events(time_limit=1.0)
                # ch._connection.sleep(3)

        # ex) WORKER_LIFE_LEFT_IN_MINUTE = 20min and TASK_KILL_TIMEOUT_MINUTE = 10min
        # timeout will be set as 10min
        # TASK_KILL_TIMEOUT_MINUTE is a extra housekeeping time after explicitly
        # terminate a task.
        logger.debug("worker life: %d", WORKER_LIFE_LEFT_IN_MINUTE.value)
        time_out_in_minute = int(WORKER_LIFE_LEFT_IN_MINUTE.value - TASK_KILL_TIMEOUT_MINUTE)
        logger.info("Timeout in minute: %d", time_out_in_minute)

    proc_return_code = -1

    logger.debug("Start subprocess to run a task.")
    try:
        p = subprocess.Popen(user_task_cmd.split(),
                             env=os.environ,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
    except MemoryError:
        logger.exception("Exception: Out of memory, %s", user_task_cmd)
        proc_return_code = DONE_FLAGS["failed with out-of-mem"]
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
            proc_return_code = DONE_FLAGS["failed with timeout"]
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
        if proc_return_code == DONE_FLAGS["failed with timeout"]:
            logger.info("User task timeout. Not enough time left for the worker.")
            return_msg["done_flag"] = str(DONE_FLAGS["failed with timeout"])
            return_msg["ret_msg"] = "User task timeout"
        elif proc_return_code == DONE_FLAGS["failed with out-of-mem"]:
            logger.info("User task out-of-mem.")
            return_msg["done_flag"] = str(DONE_FLAGS["failed with out-of-mem"])
            return_msg["ret_msg"] = "User task out-of-mem"
        elif proc_return_code == 2:  # system code
            logger.info("input file or command not found.")
            return_msg["done_flag"] = DONE_FLAGS["failed with input file or command not found"]
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
                           num_cores,
                           job_time, clone_time_rate, task_queue_name, pool_name,
                           nwpn, exch_name, worker_hb_queue):
    """
    Send heartbeats to the client
    :param interval: time interval to send heartbeats to the client
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
            #
            # logger.debug("rootpid = {} childpid = {}".format(PARENT_PROCESS_ID,
            #                                                  USER_PROC_PROC_ID.value))
            if USER_PROC_PROC_ID.value == 0:
                child_pid = PARENT_PROCESS_ID
            else:
                child_pid = int(USER_PROC_PROC_ID.value)

            # Collect pids from process tree
            root_pid_num = get_pid_tree(PARENT_PROCESS_ID)
            # logger.debug("get_pid_tree(root_proc_id={}) = {}".format(PARENT_PROCESS_ID,
            #                                                          root_pid_num))

            if PARENT_PROCESS_ID != child_pid:
                if child_pid == -9:  # if it's terminated by "kill"
                    child_pid = PARENT_PROCESS_ID
                pid_list_child = get_pid_tree(child_pid)
                # logger.debug("get_pid_tree(child_pid={}) = {}".format(child_pid,
                #                                                       pid_list_child))
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
            num_workers_on_node = get_num_workers_on_node(CONFIG)

            # Check if there is any task id in the ipc pipe
            if PIPE_TASK_ID_RECV.poll():
                task_id = PIPE_TASK_ID_RECV.recv()

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
                    #                             % (slurm_job_id), log=logger, stdoutPrint=False)
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
                    WORKER_LIFE_LEFT_IN_MINUTE.value = int(divmod(delta.total_seconds(), 60)[0])
            except Exception as e:
                logger.critical("Something wrong in computing remaining wall clock time: %s", e)
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

            # logger.debug("Send HB to the client: {}".format(msg_dict_to_send))

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


# -----------------------------------------------------------------------
def recv_hb_from_client_proc2(task_queue, exch_name, cl_hb_q_postfix):
    """

    :param task_queue:
    :param exch_name:
    :param cl_hb_q_postfix:
    :return:
    """
    rmq_conn = RmqConnectionHB(config=CONFIG)
    conn = rmq_conn.open()
    ch = conn.channel()

    # Declare exchange
    ch.exchange_declare(exchange=exch_name,
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

    def callback(ch, method, properties, body):
        msg_unzipped = json.loads(zloads(body))
        if msg_unzipped["task_type"] == TASK_TYPE["term"] and msg_unzipped["task_queue"] == task_queue:
            raise OSError(2, 'Worker termination request received.')
        elif msg_unzipped["task_type"] == TASK_TYPE["hb"]:
            pass

    ch.basic_consume(queue=client_hb_queue_name,
                     on_message_callback=callback,
                     auto_ack=True)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        ch.stop_consuming()
        raise


# -------------------------------------------------------------------------------
def recv_task_kill_request_proc():
    """
    Wait for task termination request from JTM

    :param exch_name:
    :param queue_name:
    :return:
    """

    rmq_conn = RmqConnectionHB(config=CONFIG)
    conn = rmq_conn.open()
    ch = conn.channel()
    exch_name = JTM_TASK_KILL_EXCH
    queue_name = JTM_TASK_KILL_Q
    routing_key = UNIQ_WORKER_ID
    logger.debug("JTM_TASK_KILL_EXCH: %s" % JTM_TASK_KILL_EXCH)
    logger.debug("JTM_TASK_KILL_Q: %s" % JTM_TASK_KILL_Q)
    logger.debug("UNIQ_WORKER_ID: %s" % UNIQ_WORKER_ID)

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
    rmq_conn = RmqConnectionHB(config=CONFIG)
    conn = rmq_conn.open()
    ch = conn.channel()

    ch.exchange_declare(exchange=exch_name,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)
    ch.queue_declare(queue=queue_name,
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
                so, se, ec = run_sh_command(jtm_worker_cmd, log=logger)
                run_sh_command(jtm_worker_cmd + " --dry-run", log=logger)

            ch.basic_ack(delivery_tag=method.delivery_tag)

        else:
            # ch.basic_nack(delivery_tag=method.delivery_tag)
            # If UNIQ_WORKER_ID is not for me, reject and requeue it
            logger.info("NACK sent to the broker")
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
    while True:
        logger.debug("Total Number of processes of the worker = %d" % (len(pid_list)+2))
        if len(pid_list) != NUM_WORKER_PROCS - 2:
            raise OSError(2, 'Number of processes is wrong')
        time.sleep(NUM_PROCS_CHECK_INTERVAL)


# -------------------------------------------------------------------------------
def ack_message(ch, reply_to, correlation_id, delivery_tag, response):
    """
    :param ch:
    :param delivery_tag:
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
        lost_connection_response["done_flag"] = str(DONE_FLAGS["failed with lost connection"])
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

    ######################################
    # OLD
    run_user_task(msg_unzipped, result_dict, ch)

    # NEW
    # Note: to fix connection lost
    # ex) 2019-06-10 12:13:37,917 | jtm-worker | do_work | CRITICAL : Something wrong in do_work():
    # Stream connection lost: error(104, 'Connection reset by peer')
    #
    # https://stackoverflow.com/questions/14572020/handling-long-running-tasks-in-pika-rabbitmq
    #
    # thread = threading.Thread(target=run_user_task,
    #                           args=(msg_unzipped, result_dict, ch))
    # thread.start()
    #
    # Note: if use this is_alive() with threaded ack method,
    # ref) https://github.com/pika/pika/blob/master/examples/basic_consumer_threaded.py
    # ==> psutil._exceptions.NoSuchProcess: psutil.NoSuchProcess process no longer exists (pid=97128)
    # So if threaded ack method is used, just call run_user_task()
    #
    # while thread.is_alive():  # Loop while the thread is processing
    #     ch._connection.sleep(1.0)
    ######################################

    json_data = json.dumps(result_dict)
    logger.debug("Reply msg with result: %s" % str(json_data))
    logger.debug(json_data)
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
           custom_job_log_dir_name: str, pool_name_param: str, dry_run: bool,
           slurm_job_id_param: int, worker_type_param: str, cluster_name_param: str,
           worker_clone_time_rate_param: float, num_workers_per_node_param: int,
           worker_id_param: str, charging_account_param: str,
           num_nodes_to_request_param: int, num_cores_to_request_param: int,
           constraint_param: str, mem_per_node_to_request_param: str,
           mem_per_cpu_to_request_param: str,
           qos_param: str, job_time_to_request_param: str) -> int:

    global CONFIG
    CONFIG = ctx.obj['config']
    debug = ctx.obj['debug']
    # config file has precedence
    config_debug = CONFIG.configparser.getboolean("SITE", "debug")
    if config_debug:
        debug = config_debug
    global DEBUG
    DEBUG = debug
    global WORKER_TYPE
    WORKER_TYPE = CONFIG.constants.WORKER_TYPE
    global HB_MSG
    HB_MSG = CONFIG.constants.HB_MSG
    global VERSION
    VERSION = CONFIG.constants.VERSION
    global COMPUTE_RESOURCES
    COMPUTE_RESOURCES = CONFIG.constants.COMPUTE_RESOURCES
    global TASK_TYPE
    TASK_TYPE = CONFIG.constants.TASK_TYPE
    global DONE_FLAGS
    DONE_FLAGS = CONFIG.constants.DONE_FLAGS
    global NUM_WORKER_PROCS
    NUM_WORKER_PROCS = CONFIG.constants.NUM_WORKER_PROCS
    global TASK_KILL_TIMEOUT_MINUTE
    TASK_KILL_TIMEOUT_MINUTE = CONFIG.constants.TASK_KILL_TIMEOUT_MINUTE

    global CNAME
    CNAME = CONFIG.configparser.get("SITE", "instance_name")
    global JTM_HOST_NAME
    JTM_HOST_NAME = CONFIG.configparser.get("SITE", "jtm_host_name")
    global JTM_INNER_REQUEST_Q
    JTM_INNER_REQUEST_Q = CONFIG.configparser.get("JTM", "jtm_inner_request_q")
    global CTR
    CTR = CONFIG.configparser.getfloat("JTM", "clone_time_rate")
    global JTM_INNER_MAIN_EXCH
    JTM_INNER_MAIN_EXCH = CONFIG.configparser.get("JTM", "jtm_inner_main_exch")
    global JTM_CLIENT_HB_EXCH
    JTM_CLIENT_HB_EXCH = CONFIG.configparser.get("JTM", "jtm_client_hb_exch")
    global JTM_WORKER_HB_EXCH
    JTM_WORKER_HB_EXCH = CONFIG.configparser.get("JTM", "jtm_worker_hb_exch")
    global CLIENT_HB_Q_POSTFIX
    CLIENT_HB_Q_POSTFIX = CONFIG.configparser.get("JTM", "client_hb_q_postfix")
    global WORKER_HB_Q_POSTFIX
    WORKER_HB_Q_POSTFIX = CONFIG.configparser.get("JTM", "worker_hb_q_postfix")
    global JTM_TASK_KILL_EXCH
    JTM_TASK_KILL_EXCH = CONFIG.configparser.get("JTM", "jtm_task_kill_exch")
    global JTM_TASK_KILL_Q
    JTM_TASK_KILL_Q = CONFIG.configparser.get("JTM", "jtm_task_kill_q")
    global JTM_WORKER_POISON_EXCH
    JTM_WORKER_POISON_EXCH = CONFIG.configparser.get("JTM", "jtm_worker_poison_exch")
    global JTM_WORKER_POISON_Q
    JTM_WORKER_POISON_Q = CONFIG.configparser.get("JTM", "jtm_worker_poison_q")
    global NUM_PROCS_CHECK_INTERVAL
    NUM_PROCS_CHECK_INTERVAL = CONFIG.configparser.getfloat("JTM", "num_procs_check_interval")
    global ENV_ACTIVATION
    ENV_ACTIVATION = CONFIG.configparser.get("JTM", "env_activation")
    WORKER_CONFIG_FILE = CONFIG.configparser.get("JTM", "worker_config_file")

    RMQ_HOST = CONFIG.configparser.get("RMQ", "host")
    RMQ_PORT = CONFIG.configparser.get("RMQ", "port")
    USER_NAME = CONFIG.configparser.get("SITE", "user_name")
    PRODUCTION = False
    if CONFIG.configparser.get("JTM", "run_mode") == "prod":
        PRODUCTION = True
    JOBTIME = CONFIG.configparser.get("SLURM", "jobtime")
    CONSTRAINT = CONFIG.configparser.get("SLURM", "constraint")
    CHARGE_ACCNT = CONFIG.configparser.get("SLURM", "charge_accnt")
    QOS = CONFIG.configparser.get("SLURM", "qos")
    PARTITION = CONFIG.configparser.get("SLURM", "partition")
    MEMPERCPU = CONFIG.configparser.get("SLURM", "mempercpu")
    MEMPERNODE = CONFIG.configparser.get("SLURM", "mempernode")
    NWORKERS = CONFIG.configparser.getint("JTM", "num_workers_per_node")
    NCPUS = CONFIG.configparser.getint("SLURM", "ncpus")

    global FILE_CHECK_INTERVAL
    FILE_CHECK_INTERVAL = CONFIG.configparser.getfloat("JTM", "file_check_interval")
    global FILE_CHECKING_MAX_TRIAL
    FILE_CHECKING_MAX_TRIAL = CONFIG.configparser.getint("JTM", "file_checking_max_trial")
    global FILE_CHECK_INT_INC
    FILE_CHECK_INT_INC = CONFIG.configparser.getfloat("JTM", "file_check_int_inc")

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

    print("JTM Worker, version: {}".format(VERSION))

    # Set uniq worker id if worker id is provided in the params
    if worker_id_param:
        global UNIQ_WORKER_ID
        UNIQ_WORKER_ID = worker_id_param

    # Logger setting
    log_level = "info"
    if DEBUG:
        log_level = "debug"

    setup_custom_logger(log_level, log_dir_name,
                        1, 1,
                        worker_id=UNIQ_WORKER_ID)

    logger.info("\n*****************\nDebug mode is %s\n*****************"
                % ("ON" if DEBUG else "OFF"))

    hearbeat_interval = CONFIG.configparser.getfloat("JTM", "worker_hb_send_interval")

    logger.info("Set jtm log file location to %s", log_dir_name)
    logger.info("Set jtm job file location to %s", job_script_dir_name)
    logger.info("RabbitMQ broker: %s", RMQ_HOST)
    logger.info("RabbitMQ port: %s", RMQ_PORT)
    logger.info("Pika version: %s", pika.__version__)
    logger.info("JTM user name: %s", USER_NAME)
    logger.info("Unique worker ID: %s", UNIQ_WORKER_ID)
    logger.info("\n*****************\nRun mode is %s\n*****************"
                % ("PROD" if PRODUCTION else "DEV"))
    logger.info("env activation: %s", ENV_ACTIVATION)
    logger.info("JTM config file: %s" % (CONFIG.config_file))

    # Slurm config
    num_nodes_to_request = 0
    if num_nodes_to_request_param:
        num_nodes_to_request = num_nodes_to_request_param
        # Todo
        #  Cori and JGI Cloud are exclusive allocation. So this is not needed.
        # assert mem_per_node_to_request_param is not None, "-N needs --mem-per-cpu (-mc) setting."

    # 11.13.2018 decided to remove all default values from argparse
    num_workers_per_node = num_workers_per_node_param if num_workers_per_node_param else NWORKERS
    assert num_workers_per_node > 0
    mem_per_cpu_to_request = mem_per_cpu_to_request_param if mem_per_cpu_to_request_param else MEMPERCPU
    mem_per_node_to_request = mem_per_node_to_request_param if mem_per_node_to_request_param else MEMPERNODE
    assert mem_per_cpu_to_request
    assert mem_per_node_to_request
    num_cpus_to_request = num_cores_to_request_param if num_cores_to_request_param else NCPUS
    assert num_cpus_to_request

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
                # Option 1
                # mem_per_node_to_request_byte = system_free_mem_bytes
                # Option 2
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

    # Start hb receive thread
    tp_name = ""
    if pool_name_param:
        tp_name = pool_name_param

    assert pool_name_param is not None, "User pool name is not set"
    inner_task_request_queue = JTM_INNER_REQUEST_Q + "." + pool_name_param

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

        worker_config = CONFIG.config_file if CONFIG else ""
        if WORKER_CONFIG_FILE:
            worker_config = WORKER_CONFIG_FILE

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

                        if mem_per_node_to_request:
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
                            batch_job_misc_params += " -wi %(worker_id)s_${i}" % \
                                                     dict(worker_id=UNIQ_WORKER_ID)

                    ###########################
                    if 1:
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
%(export_jtm_config_file)s
for i in {1..%(num_workers_per_node)d}
do
    echo "jobid: $SLURM_JOB_ID"
    jtm %(set_jtm_config_file)s %(debug)s worker --slurm_job_id $SLURM_JOB_ID \
-cl cori \
-wt %(worker_type)s \
-t %(wall_time)s \
--clone_time_rate %(clone_time_rate)f %(task_queue)s \
--num_worker_per_node %(num_workers_per_node)d \
-C %(constraint)s \
-m %(mem)s \
%(other_params)s &
    sleep 1
done
wait
""" % \
                                                dict(debug="--debug" if DEBUG else "",
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
                                                     mem=mem_per_node_to_request,
                                                     job_name=job_name,
                                                     exclusive=excl_param,
                                                     export_jtm_config_file="export JTM_CONFIG_FILE=%s"
                                                                            % worker_config,
                                                     set_jtm_config_file="--config=%s"
                                                                         % worker_config)

                elif cluster_name in ("lawrencium", "jgi_cloud", "jaws_lbl_gov", "jgi_cluster", "lbl"):

                    if worker_id_param:
                        batch_job_misc_params += " -wi %(worker_id)s_${i}" \
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
-cl %(lbl_cluster_name)s \
-wt %(worker_type)s \
-t %(wall_time)s \
--clone_time_rate %(clone_time_rate)f %(task_queue)s \
--num_worker_per_node %(num_workers_per_node)d \
-m %(mem)s \
%(other_params)s &
    sleep 1
done
wait
""" % \
                                            dict(debug="--debug" if DEBUG else "",
                                                 wall_time=job_time_to_request,
                                                 job_name=job_name,
                                                 partition_name=part_param,
                                                 qosname=qos_param,
                                                 charging_account=charge_param,
                                                 num_nodes_to_request=nnode_param,
                                                 mem_per_node_setting=mnode_param,
                                                 worker_id=UNIQ_WORKER_ID,
                                                 job_dir=job_script_dir_name,
                                                 env_activation_cmd=ENV_ACTIVATION,
                                                 num_workers_per_node=num_workers_per_node,
                                                 mem=mem_per_node_to_request,
                                                 lbl_cluster_name=cluster_name,
                                                 worker_type=THIS_WORKER_TYPE,
                                                 clone_time_rate=worker_clone_time_rate,
                                                 task_queue=tp_param,
                                                 other_params=batch_job_misc_params,
                                                 export_jtm_config_file="export JTM_CONFIG_FILE=%s"
                                                                        % worker_config,
                                                 set_jtm_config_file="--config=%s"
                                                                     % worker_config)

                jf.writelines(batch_job_script_str)

            os.chmod(batch_job_script_file, 0o775)

            if dry_run:
                print(batch_job_script_str)
                sys.exit(0)

            sbatch_cmd = "sbatch --parsable %s" % (batch_job_script_file)
            _, _, ec = run_sh_command(sbatch_cmd, log=logger)
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
    logger.debug("Main pid = {}".format(PARENT_PROCESS_ID))

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
        recv_hb_from_client_proc_hdl = mp.Process(target=send_hb_to_client_proc,
                                                  args=(hearbeat_interval,
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
        recv_hb_from_client_proc_hdl.start()
        pid_list.append(recv_hb_from_client_proc_hdl)
    except Exception as e:
        logger.exception("send_hb_to_client_proc: {}".format(e))
        proc_clean(pid_list)
        conn_clean(conn, ch)
        sys.exit(1)

    logger.info("Start sending my heartbeat to the client in every %d sec to %s"
                % (hearbeat_interval, WORKER_HB_Q_POSTFIX))

    # Start poison receive thread
    try:
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
        recv_poison_proc_hdl.start()
        pid_list.append(recv_poison_proc_hdl)
    except Exception as e:
        logger.exception("recv_reproduce_or_die_proc: {}".format(e))
        proc_clean(pid_list)
        conn_clean(conn, ch)
        sys.exit(1)

    # Start hb send thread
    try:
        send_hb_to_client_proc_hdl = mp.Process(target=recv_hb_from_client_proc2,
                                                args=(inner_task_request_queue,
                                                      JTM_CLIENT_HB_EXCH,
                                                      CLIENT_HB_Q_POSTFIX))
        send_hb_to_client_proc_hdl.start()
        pid_list.append(send_hb_to_client_proc_hdl)
    except Exception as e:
        logger.exception("Worker termination request: {}".format(e))
        proc_clean(pid_list)
        conn_clean(conn, ch)
        sys.exit(1)

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
