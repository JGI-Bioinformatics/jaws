#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Copyright 2018-2019 (C) DOE JGI, LBL
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
    10.26.2018 1.4.8: Fixed user termination error code update (-4); Set g_userProcPid = -9 if user
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
                      custom pool name as task queue if -tp is used; Updated routine to get the
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


"""
from jtm import *

from multiprocessing import Pipe
from datetime import timedelta

# This ipc pipe is to send taskId (by on_request())to send_hb_to_client_thread()
# when a task is requested
taskIdSendP, taskIdRecvP = Pipe()

print("JTM Worker, version: {}".format(VERSION))  # VERSION <- Config.py

# -------------------------------------------------------------------------------
# Globals
# -------------------------------------------------------------------------------
g_recvHbFromClientProc = None
g_sendHbToClientProc = None
g_recvPoisonProc = None
g_taskKillProc = None

g_userProcPid = multiprocessing.Value("i", 0)
g_allowZeroOutputFile = False  # False ==> if size(output file)==0, return error code
g_uniqWorkerId = str(shortuuid.uuid())
g_bClientAlive = False
g_workerStartTime = datetime.datetime.now()


# -------------------------------------------------------------------------------
def run_something(msgUnzipped, msgContainer):
    """
    Run a user command in msgZipped
    :param msgUnzipped: compressed msg from client
    :return:
    """
    # Uncompress msg to get a task
    taskId = msgUnzipped["task_id"]
    userCmd = msgUnzipped["user_cmd"]
    outFiles = msgUnzipped["output_files"]
    taskType = msgUnzipped["task_type"]
    # cromwellJid = msgUnzipped["cromwell_jid"]

    # msgContainer = {}
    msgContainer["task_id"] = taskId
    msgContainer["user_cmd"] = userCmd
    msgContainer["task_type"] = taskType

    if "stdout" in msgUnzipped:
        userCmd + " > %s" % (msgUnzipped["stdout"])
    if "stderr" in msgUnzipped:
        userCmd + " 2>%s" % (msgUnzipped["stderr"])

    # Run the task
    logger.info("Running task ID %d...", taskId)
    # logger.info("Running task/task ID, %s/%d...", cromwellJid, taskId)
    try:
        p = subprocess.Popen(userCmd,
                             shell=True,
                             env=os.environ,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)

        # Set g_userProcPid = forked child process id (this value will be sent
        # to send_hb_to_client function)
        g_userProcPid.value = p.pid
        logger.debug("g_userProcPid.value %d" % g_userProcPid.value)
    except Exception as detail:
        logger.exception("Exception: Failed to run user command, %s", userCmd)
        logger.exception("Detail: %s", str(detail))
        raise

    logger.debug("User command: %s", userCmd)
    stdoutValue = p.communicate()[0]
    p.wait()
    p.poll()

    # Prepare result to send back
    msgContainer["out_files"] = outFiles
    msgContainer["worker_id"] = g_uniqWorkerId
    msgContainer["host_name"] = socket.gethostname()

    if g_userProcPid.value == -9:
        msgContainer["done_flag"] = DONE_FLAGS["failed with user termination"]
        msgContainer["ret_msg"] = "Task cancelled."
    else:
        if p.returncode == 2:
            logger.info("input file not found.")
            msgContainer["done_flag"] = "1"
            msgContainer["ret_msg"] = "Input file not found."

        elif p.returncode == 0:
            logger.info("Task# %s completed!" % (str(taskId)))

            # Output file checking
            if len(outFiles):
                ofs = outFiles.split(",")
                logger.debug("Number of output files = %d.", len(ofs))
                outFileList = []

                for i in range(len(ofs)):
                    outFileList.append(ofs[i])

                ret, fSize = check_output(outFileList,
                                          FILE_CHECK_INTERVAL,
                                          FILE_CHECKING_MAX_TRIAL,
                                          FILE_CHECK_INT_INC)

                if not ret:
                    retMsg = "Failed to check output file(s): %s, file size = %s." % (ofs, fSize)
                    logger.critical(retMsg)
                    msgContainer["done_flag"] = DONE_FLAGS["failed to check output file(s)"]
                    msgContainer["ret_msg"] = retMsg
                else:
                    msgContainer["done_flag"] = DONE_FLAGS["success with correct output file(s)"]
                    msgContainer["ret_msg"] = "Output file checking is OK."
            else:
                msgContainer["done_flag"] = DONE_FLAGS["success"]
                msgContainer["ret_msg"] = "No file(s) to check."
        else:
            logger.critical("Failed to execute a task, %s. Non-zero exit code. stdout = %s."
                            % (userCmd, stdoutValue))
            msgContainer["done_flag"] = DONE_FLAGS["failed to run user command"]
            msgContainer["ret_msg"] = stdoutValue

    # jsonData = json.dumps(msgContainer)
    # logger.debug("Result reply: %s" % str(jsonData))
    # msgZippedToSend = zdumps(jsonData)
    # return msgZippedToSend


# -------------------------------------------------------------------------------
def ack_message(ch, deliveryTag):
    """

    :param ch:
    :param delivery_tag:
    :return:
    """

    # Note that `ch` must be the same pika channel instance via which
    #  the message being ACKed was retrieved (AMQP protocol constraint).

    if ch.is_open:
        ch.basic_ack(deliveryTag)
    else:
        # Channel is already closed, so we can't ACK this message;
        # log and/or do something that makes sense for your app in this case.
        pass


# -------------------------------------------------------------------------------
def send_result_back_to_manger(response, replyTo, corrId):
    """

    :param response:
    :param replyTo:
    :param corrId:
    :return:
    """
    try:
        rabbitConnection = RmqConnectionHB()
        conn = rabbitConnection.open()
        ch = conn.channel()
        ch.queue_declare(queue=replyTo,
                         durable=True,
                         exclusive=False,
                         auto_delete=True)
        assert replyTo.endswith(CNAME)
        ch.basic_publish(exchange=JTM_INNER_MAIN_EXCH,
                         routing_key=replyTo,  # use the queue which the client created
                         properties=pika.BasicProperties(
                             delivery_mode=2,  # make message persistent
                             correlation_id=corrId),
                         body=response)
    except Exception as e:
        logger.critical("Something wrong in send_result_back_to_manger(): %s", e)
        raise

    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def do_task_processing(conn, ch, deliveryTag, replyTo, corrId, body):
    """

    :param conn:
    :param ch:
    :param deliveryTag:
    :param replyTo:
    :param corrId:
    :param body:
    :return:
    """
    # thread_id = threading.get_ident()
    # logger.info('Thread id: %s Delivery tag: %s Message body: %s', thread_id, deliveryTag, body)

    msgUnzipped = json.loads(zloads(body))
    taskId = msgUnzipped["task_id"]
    taskIdSendP.send(taskId)
    if taskId == -9:
        logger.info("Received TERM signal. Terminate myself.")
        #
        # Ref) http://gavinroy.com/deeper-down-the-rabbit-hole-of-message-redeli
        # Return back the message (the TERM message) in the queue
        # so that the other workers can get it.
        #
        # ch.basic_reject(delivery_tag=method.delivery_tag, requeue=True)
        ch.basic_reject(delivery_tag=deliveryTag, requeue=True)
        ch.stop_consuming()
        ch.close()

        if g_recvHbFromClientProc: g_recvHbFromClientProc.terminate()
        if g_sendHbToClientProc: g_sendHbToClientProc.terminate()
        if g_recvPoisonProc: g_recvPoisonProc.terminate()
        if g_taskKillProc: g_taskKillProc.terminate()

        sys.exit(0)

    logger.info("Received a task, %r" % (msgUnzipped,))
    logger.debug("Return queue = %s", replyTo)

    # response = run_something(msgUnzipped)
    resDict = {}
    run_something(msgUnzipped, resDict)
    # thread = threading.Thread(target=run_something, args=(msgUnzipped, resDict))
    # thread.start()
    # while thread.is_alive():  # Loop while the thread is processing
    #     ch._connection.sleep(1.0)
    jsonData = json.dumps(resDict)
    logger.debug("Result reply: %s" % str(jsonData))
    response = zdumps(jsonData)

    try:
        logger.debug("Sent the result back to the client via '%s' queue.", replyTo)
        send_result_back_to_manger(response, replyTo, corrId)

        # Note: After sending ack, the message will be deleted from RabbitMQ
        #  If this worker crashes while running a user command, this task will
        #  be sent to other workers available
        # ch.basic_ack(delivery_tag=method.delivery_tag)
        cb = functools.partial(ack_message, ch, deliveryTag)
        conn.add_callback_threadsafe(cb)

    except Exception as e:
        logger.critical("Something wrong in on_request(): %s", e)

    # Send taskid=0 to send_hb_to_client_thread() b/c the requested task is completed
    taskIdSendP.send(0)

    # Reset child pid
    g_userProcPid.value = 0


# -------------------------------------------------------------------------------
# def on_task_request(ch, method_frame, _header_frame, body, args):
#     """
#     :param ch:
#     :param method_frame:
#     :param _header_frame:
#     :param body:
#     :param args:
#     :return:
#     """
#     (conn, thrds) = args
#     # print method_frame, _header_frame
#     delivery_tag = method_frame.delivery_tag
#     reply_to = _header_frame.reply_to
#     correlation_id = _header_frame.correlation_id
#     t = threading.Thread(target=do_task_processing,
#                          args=(conn,
#                                ch,
#                                delivery_tag,
#                                reply_to,
#                                correlation_id,
#                                body))
#     t.start()
#     thrds.append(t)


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
    msgUnzipped = json.loads(zloads(body))
    taskId = msgUnzipped["task_id"]

    # Get the length of the queue
    # taskQueueLen = ch.queue_declare(queue=method.routing_key,
    #                                 durable=True,
    #                                 exclusive=False,
    #                                 auto_delete=True).method.message_count
    # logger.debug("#tasks queued = %d", taskQueueLen)

    # Send taskid to send_hb_to_client_thread() so that the taskId can be shown in the hb messages
    taskIdSendP.send(taskId)

    # If client sends "terminate worker"
    if taskId == -9:
        logger.info("Received TERM signal. Terminate myself.")
        #
        # Ref) http://gavinroy.com/deeper-down-the-rabbit-hole-of-message-redeli
        # Return back the message (the TERM message) in the queue
        # so that the other workers can get it.
        #
        ch.basic_reject(delivery_tag=method.delivery_tag, requeue=True)
        ch.stop_consuming()
        ch.close()

        if g_recvHbFromClientProc: g_recvHbFromClientProc.terminate()
        if g_sendHbToClientProc: g_sendHbToClientProc.terminate()
        if g_recvPoisonProc: g_recvPoisonProc.terminate()
        if g_taskKillProc: g_taskKillProc.terminate()

        sys.exit(0)

    logger.info("Received a task, %r" % (msgUnzipped,))
    logger.debug("Return queue = %s", props.reply_to)

    ######################################
    # response = run_something(msgUnzipped)

    # Note: to fix connection lost
    #  2019-06-10 12:13:37,917 | jtm-worker | on_request | CRITICAL : Something wrong in on_request(): Stream connection lost: error(104, 'Connection reset by peer')
    resDict = {}
    thread = threading.Thread(target=run_something, args=(msgUnzipped, resDict))
    thread.start()
    while thread.is_alive():  # Loop while the thread is processing
        ch._connection.sleep(1.0)
    jsonData = json.dumps(resDict)
    logger.debug("Result reply: %s" % str(jsonData))
    response = zdumps(jsonData)
    ######################################

    try:
        logger.debug("Sent the result back to the client via '%s' queue.", props.reply_to)

        # Note: to keep result messages even the client is not alive.
        ch.queue_declare(queue=props.reply_to,
                         durable=True,
                         exclusive=False,
                         auto_delete=True)
                         # auto_delete=False)

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
    taskIdSendP.send(0)

    # Reset child pid
    g_userProcPid.value = 0


# -------------------------------------------------------------------------------
def check_output(outFiles, waitSecOrig=3, maxTrials=3, waitSecIncrease=1.5):
    """
    Check 1) existence, 2) size>0 for each file in outFiles
    :param outFiles: list of absolute paths to the output files to check
    :param waitSecOrig: sleep time between output file checking before retiral
    :param maxTrials: max trial for checking
    :param waitSecIncrease: wait time increase for retrial
    :return:
    """
    fSize = 0
    fOk = False
    trial = 1

    for ofile in outFiles:
        logger.info("Output file check: %s", ofile)
        logger.debug("Output file check: %s", os.path.expandvars(ofile))
        ofile = os.path.expandvars(ofile)

        fSize = 0
        fOk = False

        while trial < maxTrials:
            logger.info("Output file checking. Trial# = %d", trial)

            # First, check file existence
            # os.path.exists returns if it is a valid path(check for directory or file, both)
            # and os.path.isfile(checks for only file, not directory) returns if it is a file
            bFileExist = os.path.isfile(ofile)

            # If exist, check file size
            if bFileExist:
                fSize = os.path.getsize(ofile)
                if fSize == 0:
                    logger.warning("File, %s is zero size.", ofile)

                # If "-z" option is used, allow zero sized output file
                if g_allowZeroOutputFile and fSize == 0:
                    fSize = 1

                if bFileExist and fSize > 0:
                    fOk = True
                    logger.info ("Output file '%s' is OK.", ofile)
                    break
            else:
                logger.info("Outout file not found.")

            # Wait for initial wait time
            time.sleep(waitSecOrig)

            # Increase the wait time
            waitSecOrig *= waitSecIncrease
            trial += 1

    return fOk, fSize


# -------------------------------------------------------------------------------
def send_hb_to_client_thread(rootPid, interval, qChildPid, workerId, aPipe, slurmJobId, workerType,
                             mpn, mpc, ncores, jobtime, ctr, taskQ, poolName, nwpn):
    """
    Send heartbeats to the client
    :param rootPid: rootPid of this worker
    :param interval: time interval to send heartbeats to the client
    :param qChildPid: child pid for user command
    :param workerId: worker id
    :param aPipe: Pythin IPC pipe for getting task id
    :param slurmJobId: SLURM job id
    :param workerType: worker type [manual, static, dynamic]
    :param mpn: memory request per node
    :param mpc: memory request per core
    :param ncores: number of cores
    :param jobtime: wallclocktime
    :param ctr: clone time rate
    :param taskQ: task queue name
    :param poolName: pool name
    :return:
    """
    # Remote broker (mq.nersc.gov)
    # rabbitConnection = RmqConnection()
    # with heartbeat_interval=0
    # ref) http://stackoverflow.com/questions/14572020/handling-long-running-tasks-in-pika-rabbitmq,
    # http://stackoverflow.com/questions/34721178/pika-blockingconnection-rabbitmq-connection-closed
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()

    exchName = JTM_WORKER_HB_EXCH
    workerHbQ = WORKER_HB_Q_POSTFIX
    hostName = socket.gethostname()

    ch.exchange_declare(exchange=exchName,
                        exchange_type="direct",
                        passive=False,
                        durable=False,
                        auto_delete=False)

    taskId = 0
    pidListMerged = []

    while 1:
        try:
            # TODO: make it optional to send the resource info on hb message
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
            logger.debug("rootpid = {} childpid = {}".format(rootPid, qChildPid))
            rootPid = int(rootPid)
            if qChildPid.value == 0:
                childPid = rootPid
            else:
                childPid = int(qChildPid.value)

            # Collect pids from process tree
            pidListRoot = get_pid_tree(rootPid)
            logger.debug("get_pid_tree(rootPid={}) = {}".format(rootPid, pidListRoot))

            if rootPid != childPid:
                if childPid == -9:  # if it's terminated by "kill"
                    childPid = rootPid
                pidListChild = get_pid_tree(childPid)
                logger.debug("get_pid_tree(childPid={}) = {}".format(childPid, pidListChild))
                if len(pidListChild) > 0:
                    pidListMerged = pidListRoot + pidListChild[1:]
                else:
                    # NOTE: be careful on this resetting! Might lose child pid
                    # qChildPid.value = 0
                    pidListMerged = pidListRoot
            else:
                pidListMerged = pidListRoot

            vmemUsageList = []
            rmemUsageList = []
            try:
                vmemUsageList.extend([get_virtual_memory_usage(pid, 0.0, False) for pid in pidListMerged])
            except ValueError:
                logger.exception("ValueError: Failed to collect VM memory usage.")
                pass
            except UnboundLocalError:
                logger.exception("UnboundLocalError: No entry in process id list.")
                pass

            try:
                rmemUsageList.extend([get_resident_memory_usage(pid, 0.0, False) for pid in pidListMerged])
            except ValueError:
                logger.exception("ValueError: Failed to collect RES memory usage.")
                pass
            except UnboundLocalError:
                logger.exception("UnboundLocalError: No entry in process id list.")
                pass

            # Collect cpu_usages for all pids in the tree and get max()
            cpuLoadList = [get_cpu_load(pid) for pid in pidListMerged]
            cpuLoad = max(cpuLoadList) if len(cpuLoadList) > 0 else 0.0

            # Collect mem_usages for all pids in the tree and get sum()
            rmemUsage = "%.1f" % sum(rmemUsageList)
            vmemUsage = "%.1f" % sum(vmemUsageList)

            # Only get the run time of childPid
            if sys.platform.lower() == "darwin":
                # TODO: Add a method to get etime on Mac OS
                runTime = ""
            else:
                runTime = get_runtime(childPid)

            if cpuLoad == "":
                cpuLoad = 0.0

            if runTime == "":
                runTime = 0

            # Get % mem used per node
            # This is for node-based scheduling
            percUsedMem = "%.1f" % get_total_mem_usage_per_node()
            numWorkersOnThisNode = get_num_workers_on_node()

            # Check if there is any task id in the ipc pipe
            if aPipe.poll():
                taskId = aPipe.recv()

            today = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            endDatetime = None
            lifeInMinute = 0

            try:
                # Note: this only works with SLURM. For other HPCs/Clouds, methods to get the
                #  remaining life are needed! ==> resovled!
                if slurmJobId:
                    # OLD
                    # fixme; slurm overload
                    # so, se, ec = run_sh_command("scontrol show jobid %d"
                    #                             % (slurmJobId), live=True, log=logger, stdoutPrint=False)
                    # pat = re.compile(r"""([^\s=]+)=\s*((?:[^\s=]+(?:\s|$))*)""")
                    # # Convert a list of x=y to dict
                    # entries = dict((k, v.split()) for k, v in pat.findall(so))
                    # delta = datetime.datetime.now() - datetime.datetime.strptime(entries["EndTime"][0], "%Y-%m-%dT%H:%M:%S")
                    # lifeInMinute = -divmod(delta.total_seconds(), 60)[0]
                    # endDatetime = entries["EndTime"][0]

                    # NEW
                    # jobtime2 = datetime.datetime.strptime(jobtime, '%H:%M:%S')
                    # endDatetime = g_workerStartTime + timedelta(seconds=jobtime2.second+jobtime2.minute*60+jobtime2.hour*3600)
                    # 24:00:00 --> seconds
                    jobtime2Seconds = int(jobtime.split(':')[0]) * 3600 + int(jobtime.split(':')[1]) * 60 + int(jobtime.split(':')[2])
                    endDatetime = g_workerStartTime + timedelta(seconds=jobtime2Seconds)
                    delta = endDatetime - datetime.datetime.now()
                    lifeInMinute = divmod(delta.total_seconds(), 60)[0]
            except Exception as e:
                logger.critical("Something wrong in computing remaining wall clock time: %s", e)
                # raise
                sys.exit(1)

            # Get the total number of tasks in the task queue for this pool
            taskQueueLen = ch.queue_declare(queue=taskQ,
                                            durable=True,
                                            exclusive=False,
                                            auto_delete=True).method.message_count
            # logger.debug("#tasks queued=%d in %s" % (taskQueueLen, taskQ))

            # To decrease the hb message size, changed to int key value
            msgDict = {HB_MSG["child_pid"]: childPid,
                       HB_MSG["clone_time_rate"]: ctr,
                       HB_MSG["cpu_load"]: cpuLoad,
                       # HB_MSG["end_date"]: endDatetime.strftime("%Y-%m-%dT%H:%M:%S") if endDatetime else today,
                       HB_MSG["end_date"]: today,  # Note: for dynamic worker endDate update
                       HB_MSG["host_name"]: hostName,
                       HB_MSG["ip_address"]: socket.gethostbyname(hostName),
                       HB_MSG["job_time"]: jobtime if WORKER_TYPE[workerType] > 0 else None,
                       HB_MSG["jtm_host_name"]: JTM_HOST_NAME,
                       HB_MSG["life_left"]: lifeInMinute,
                       HB_MSG["mem_per_core"]: mpc if WORKER_TYPE[workerType] > 0 else "",
                       HB_MSG["mem_per_node"]: mpn if WORKER_TYPE[workerType] > 0 else "",
                       HB_MSG["num_cores"]: ncores if WORKER_TYPE[workerType] > 0 else multiprocessing.cpu_count(),
                       HB_MSG["num_tasks"]: taskQueueLen,
                       HB_MSG["num_workers_on_node"]: numWorkersOnThisNode,
                       HB_MSG["perc_mem_used"]: percUsedMem,
                       # HB_MSG["pool_name"]: taskQ if taskQ[0] != '_' else "",  # TODO: This might cause error
                       HB_MSG["pool_name"]: poolName,
                       HB_MSG["ret_msg"]: "hb",
                       HB_MSG["rmem_usage"]: rmemUsage,
                       HB_MSG["root_pid"]: rootPid,
                       HB_MSG["run_time"]: runTime,
                       HB_MSG["slurm_jobid"]: slurmJobId,
                       HB_MSG["task_id"]: taskId,
                       HB_MSG["vmem_usage"]: vmemUsage,
                       HB_MSG["worker_id"]: workerId,
                       HB_MSG["worker_type"]: WORKER_TYPE[workerType],
                       HB_MSG["nwpn"]: nwpn}

            msgZipped = zdumps(json.dumps(msgDict))

            assert workerHbQ.endswith(CNAME)
            ch.basic_publish(exchange=exchName,
                             routing_key=workerHbQ,
                             body=msgZipped)

            logger.debug("Send HB to the client: {}".format(msgDict))

        except Exception as e:
            logger.critical("Something wrong with send_hb_to_client(): %s", e)
            os._exit(1)

        # TODO: Need to start with shorted interval like (0.5-1 sec) for about 30 sec
        #  from the beginning and then use the user specified interval
        # time.sleep(float(interval))
        conn.process_data_events(time_limit=float(interval))

    # unreachable
    ch.close()
    conn.close()


# -----------------------------------------------------------------------
def recv_hb_from_client_thread(tqueue, timeOut, ppid, workerId):
    """
    Receive heartbeat from the client
    :param tqueue:
    :param timeOut: time out in sec (if no heartbeat for timeOut, will kill myself)
    :param ppid: parent process id for termination
    :param workerId:
    :return:
    """
    # TODO: change to event driven by ch.basic_consume and ch.start_consuming
    #  Remote broker (mq.nersc.gov)
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()

    exchName = JTM_CLIENT_HB_EXCH
    hbInt = WORKER_HB_RECV_INTERVAL
    hname = socket.gethostname()

    # Declare exchange
    ch.exchange_declare(exchange=exchName,
                        # exchange_type="fanout",
                        exchange_type="topic",
                        durable=False,
                        auto_delete=False)

    # Declare queue
    clHbQueueName = "_jtm_worker_%s%s" % (workerId, CLIENT_HB_Q_POSTFIX)
    ch.queue_declare(queue=clHbQueueName,
                     durable=False,
                     exclusive=True,
                     auto_delete=True)

    # Bind exchange to hb queue for receiving the client hb
    ch.queue_bind(exchange=exchName,
                  queue=clHbQueueName,
                  routing_key="*." + CNAME)

    def panic():
        logger.info("No heartbeat from the client for %s sec. Terminate myself!", timeOut)
        if g_recvHbFromClientProc: g_recvHbFromClientProc.terminate()
        if g_sendHbToClientProc: g_sendHbToClientProc.terminate()
        if g_recvPoisonProc: g_recvPoisonProc.terminate()
        if g_taskKillProc: g_taskKillProc.terminate()
        ch.stop_consuming()
        os.kill(int(ppid), signal.SIGTERM)
        sys.exit(1)

    logger.info("Listen to %s queue for the client's heartbeat.", clHbQueueName)

    # Callback version
    # def on_client_hb(ch, method, props, body):
    #     logger.info("Client hb received: {} {} {}".format(method, props, body))
    #
    # ch.basic_qos(prefetch_count=1)  # Not to assign more than 1 message
    # ch.basic_consume(on_client_hb, queue=clHbQueueName, auto_ack=True)
    # ch.start_consuming()

    def stop():
        logger.info("Terminate myself: {} on {}".format(workerId, hname))
        if g_recvHbFromClientProc: g_recvHbFromClientProc.terminate()
        if g_sendHbToClientProc: g_sendHbToClientProc.terminate()
        if g_recvPoisonProc: g_recvPoisonProc.terminate()
        if g_taskKillProc: g_taskKillProc.terminate()
        ch.stop_consuming()
        os.kill(int(ppid), signal.SIGTERM)  # term jtm-worker process
        sys.exit(10)

    noAck = False

    try:
        while 1:
            if int(PIKA_VER[0]) < 1:  # v0.13.1
                methodFrame, headerFrame, body = ch.basic_get(queue=clHbQueueName, no_ack=True)
            else:  # v1.0.1 or higher
                methodFrame, headerFrame, body = ch.basic_get(queue=clHbQueueName, auto_ack=True)

            if int(timeOut) != 0 and methodFrame is None and headerFrame is None and body is None:
                if noAck is False:
                    noAck = True
                    # time.sleep(float(timeOut))
                    conn.process_data_events(time_limit=float(timeOut))
                    continue
                else:
                    panic()

            elif methodFrame is not None and headerFrame is not None and body is not None:
                noAck = False
                msgUnzipped = json.loads(zloads(body))
                if "task_type" in msgUnzipped:
                    # worker termination
                    if msgUnzipped["task_type"] == TASK_TYPE["term"] and msgUnzipped["task_queue"] == tqueue:
                        stop()
                    # hb recieved from manager
                    elif msgUnzipped["task_type"] == TASK_TYPE["hb"]:
                        pass

                # Consume all hb msgs
                for i in range(int(methodFrame.message_count)):
                    if int(PIKA_VER[0]) < 1:  # v0.13.1
                        methodFrame, headerFrame, body = ch.basic_get(queue=clHbQueueName, no_ack=True)
                    else:  # v1.0.1 or higher
                        methodFrame, headerFrame, body = ch.basic_get(queue=clHbQueueName, auto_ack=True)

                    msgUnzipped = json.loads(zloads(body))

                    if "task_type" in msgUnzipped:
                        if msgUnzipped["task_type"] == TASK_TYPE["term"] and msgUnzipped["task_queue"] == tqueue:
                            stop()
                        elif msgUnzipped["task_type"] == TASK_TYPE["hb"]:
                            # normal hb from jtm
                            pass

                logger.debug("Received heartbeats from the client.")


            # time.sleep(hbInt)
            conn.process_data_events(time_limit=float(hbInt))

    except Exception as e:
        logger.critical("Something wrong with recv_hb_from_client_thread(): %s", e)
        raise

    # unreachable
    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def recv_task_kill_request_thread(tp, wid, cl):
    """
    Wait for task termination request from JTM
    :param tp: task queue (pool) name
    :param wid: worker id
    :param cl: cluster name
    :return:
    """
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()

    exch = JTM_TASK_KILL_EXCH
    qName = JTM_TASK_KILL_Q
    rKey = str(wid)

    ch.exchange_declare(exchange=exch,
                        exchange_type="fanout",
                        durable=True,
                        auto_delete=False)
    ch.queue_declare(queue=qName,
                     # durable=False,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)
    ch.queue_bind(exchange=exch,
                  queue=qName,
                  routing_key=rKey)

    def on_kill(ch, method, props, body):
        msgUnzipped = json.loads(zloads(body))
        logger.debug("Received task termination command: %r" % msgUnzipped)

        if msgUnzipped["worker_id"] == wid:
            if msgUnzipped["child_pid"] > 0:
                logger.info("Process termination request received.")

                # This -9 is to notify run_something() that it's killed by user requests
                # Also send_hb_to_client() will check this for adjust childpid and parentpid
                # Note: this should done first to signal run_something() that the process is killed.
                g_userProcPid.value = -9

                # kill if there is child's children
                for i in get_pid_tree(msgUnzipped["child_pid"]):
                    # killCmd = "kill -9 %d" % msgUnzipped["child_pid"]
                    killCmd = "kill -9 %d" % i
                    logger.info("Executing {} for taskID, {}".format(killCmd, msgUnzipped["task_id"]))
                    so, se, ec = run_sh_command(killCmd, live=True, log=logger)
                    if ec == 0:
                        logger.info("Successfully terminate a user task process.")
                    else:
                        logger.warning("User process not found. Ignore the termination command, %s"
                                       % (killCmd))
                        # TODO: Failed to terminate a user process for some reason. How to deal
                        #  with this case?
                        # ch.basic_reject(delivery_tag=method.delivery_tag, requeue=True)

                # Kill the main child process
                # Note: can consider to use "pkill -9 -P ppid" to kill the family
                killCmd = "kill -9 %d" % msgUnzipped["child_pid"]
                logger.info("Executing {} for taskID, {}".format(killCmd, msgUnzipped["task_id"]))
                so, se, ec = run_sh_command(killCmd, live=True, log=logger)
                if ec == 0:
                    logger.info("Successfully terminate a user task process.")
                else:
                    logger.warning("User process not found. Failed to execute the command, %s" % (killCmd))

            else:
                logger.warning("No valid child process id to terminate.")

            ch.basic_ack(delivery_tag=method.delivery_tag)

        else:
            # ch.basic_nack(delivery_tag=method.delivery_tag)
            # If wid is not for me, reject and requeue it
            ch.basic_reject(delivery_tag=method.delivery_tag, requeue=True)


    # Waiting for a task kill
    ch.basic_qos(prefetch_count=1)

    if int(PIKA_VER[0]) < 1:  # v0.13.1
        ch.basic_consume(on_kill, queue=qName, no_ack=False)
    else:  # v1.0.1 or higher
        ch.basic_consume(queue=qName, on_message_callback=on_kill, auto_ack=False)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        ch.stop_consuming()


# -------------------------------------------------------------------------------
def recv_reproduce_or_die_thread(tp, wid, cl, mm, mc, nn, nc, ti, ctr, nw):
    """
    Wait for process termination request from JTM
    :param tp: task queue (pool) name
    :param wid: worker id
    :param cl: cluster name
    :param mm: memory size per node
    :param mc: memory size per core
    :param nn: number of nodes
    :param nc: number of cores
    :param ti: wallclocktime request
    :param ctr: cloning time rate (remaining runtime / total runtime)
    :param nw: number of workers per node
    :return:
    """
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()

    exch = JTM_WORKER_POISON_EXCH
    qName = JTM_WORKER_POISON_Q

    ch.exchange_declare(exchange=exch,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)
    ch.queue_declare(queue=qName,
                     # durable=False,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)
    ch.queue_bind(exchange=exch,
                  queue=qName,
                  routing_key=wid)


    def on_poison(ch, method, props, body):
        msgUnzipped = json.loads(zloads(body))

        if msgUnzipped["worker_id"] == wid:

            # If static worker cloning request
            if msgUnzipped["num_clones"] == -1:
                # Kill the worker request
                pass
            elif msgUnzipped["num_clones"] == 1:
                logger.info("Static worker cloning request received.")

                jtmWorkerCmd = "{} && jtm-worker -wt static -cl {} -t {} -ct {}".format(
                                ENV_ACTIVATION, cl, ti, ctr)

                if tp is not None:
                    jtmWorkerCmd += " -tp {}".format(tp)

                if nw != "" and int(nw) > 1:
                    jtmWorkerCmd += " -nw {}".format(nw)
                else:
                    jtmWorkerCmd += " -nw 1"

                if nn != "" and int(nn) > 0:
                    jtmWorkerCmd += " -N {} -m {} -c {}".format(nn, mm, nc)
                else:
                    jtmWorkerCmd += " -c {}".format(nc)
                    if mm not in ("", None):  # if node mem defined, set --mem
                        jtmWorkerCmd += " -m {}".format(mm)
                    else:  # if not, set mem per core
                        jtmWorkerCmd += " -mc {}".format(mc)

                logger.info("Executing {}".format(jtmWorkerCmd))
                so, se, ec = run_sh_command(jtmWorkerCmd, live=True, log=logger)
                logger.debug("{} {} {}".format(so, se, ec))
                _, _, _ = run_sh_command(jtmWorkerCmd + " --dry-run", live=True, log=logger)


            # Note: 01072019 no more auto cloning of dynamic worker
            # elif msgUnzipped["num_clones"] > 1:
            #     logger.info("Dynamic worker cloning request received.")
            #     jtmWorkerCmd = "{} && jtm-worker -wt dynamic -cl {} -t {} -ct {}".format(
            #                     ENV_ACTIVATION, cl, ti, ctr)
            #     if tp is not None:
            #         jtmWorkerCmd += " -tp {}".format(tp)
            #     if nn != "" and int(nn) > 0:
            #         jtmWorkerCmd += " -N {} -m {} -c {} -nw {}".format(nn, mm, nc, nw)
            #     else:
            #         jtmWorkerCmd += " -c {}".format(nc)
            #         if mm not in ("", None):  # if node mem defined, set --mem
            #             jtmWorkerCmd += " -m {}".format(mm)
            #         else:
            #             jtmWorkerCmd += " -mc {}".format(mc)
            #
            #     for _ in range(0, msgUnzipped["num_clones"]):
            #         logger.info("Executing {}".format(jtmWorkerCmd))
            #         _, _, _ = run_sh_command(jtmWorkerCmd, live=True, log=logger)
            #         # _, _, _ = run_sh_command(jtmWorkerCmd + " --dry-run", live=True, log=logger)

            ch.basic_ack(delivery_tag=method.delivery_tag)

        else:
            # ch.basic_nack(delivery_tag=method.delivery_tag)
            # If wid is not for me, reject and requeue it
            ch.basic_reject(delivery_tag=method.delivery_tag, requeue=True)

    # waiting for poison or cloning command
    ch.basic_qos(prefetch_count=1)
    if int(PIKA_VER[0]) < 1:  # v0.13.1
        ch.basic_consume(on_poison, queue=qName, no_ack=False)
    else:  # v1.0.1 or higher
        ch.basic_consume(queue=qName, on_message_callback=on_poison, auto_ack=False)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        ch.stop_consuming()

