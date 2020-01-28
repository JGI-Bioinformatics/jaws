#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Copyright 2018-2019 (C) DOE JGI, LBL
# Seung-Jin Sul (ssul@lbl.gov)
#
# pylint: disable=C0111,C0103,R0205
#
"""

JGI task manager


Example of task processing scenario
1. jtm-submiit sends a msg to "jgi_microservice" with "jtm_task_request_queue" tag.
2. jtm-manager listens to "jtm_task_request_queue" which is bound to "jgi_microservice"
3. When a task is published, jtm-manager takes it and sends it to a pool
   (to jgi_jtm_inner_main_exchange)
4. Workers listen to jgi_jtm_inner_request_queue which is bound to "jgi_jtm_inner_main_exchange"
5. When a task is completed, a worker sends a result msg to "jgi_jtm_inner_main_exchange" with
   "jgi_jtm_inner_result_queue" tag
6. jtm-manager listens to "jgi_jtm_inner_result_queue" queue. When a result is ready,
   takes and updates tables


Revisions:

    10.01.2015 2.6.0: Tested with pika 0.10.0

    12.03.2015 2.7.0: Tested with heartbeat_interval=0 for worker and  heartbeat_interval=60 for
                      client in BlockingConnection()
    08.22.2016 2.7.1: Checked invalid user command
    11.15.2016 2.7.2: Added "-lf" option to set jtm log file saving location
    03.03.2017 2.7.3: Used execute many to insert tasks into tasks table when "-d" option is used

    03.03.2017 3.0.0: Updated to set hearbeat_internal=0 (RmqConnectionHB(0)) so that connection
                      timeout is disalbed
    03.10.2017 3.0.1: Changed to use executemany for creating task list sqlite database
    03.16.2017 3.0.2: Added custom separator (default=':') for creating task list

    09.12.2018 3.1.0: Branched out 'jtm'; Updated for pika=0.12.0; Changed 'type' to 'exchange_type'
                      in exchange_declare()
                      Added 'jgi_jtm_task_manager' exchange for task and result messages
    09.13.2018 3.1.1: Added unique worker id

    09.13.2018 3.2.0: Added worker heartbeat exchange

    ---------------------------------------------------------------------------------------------

    09.19.2018 0.0.9: Started jgi-task-manager

    09.21.2018 1.0.0: Working version done
    09.27.2018 1.1.0: Updated message structure;

    10.03.2018 1.2.0: Added task_type; Added db utils;
    10.04.2018 1.2.1: Set jgi_jtm_client_hb_exchange durable=True and auto_delete=False so that it
                      can be maintained even with no worker; Added runs table;
    10.05.2018 1.2.2: Updated resource msg as dict
    10.05.2018 1.2.3: Updated send hb to client interval to 5sec;
                      Added comm pipe to send taskid with hb;
    10.05.2018 1.2.4: Updated to update runs table with resource data

    10.09.2018 1.3.0: Added JtmInterface class; Added jtm_submit and jtm_status
    10.10.2018 1.3.1: Changed to MySQL
    10.12.2018 1.3.2: Updated recv_hb_from_worker;
               1.3.3: Updated to make connection per each SQL;

               1.4.0: jtm_kill works
    10.15.2018 1.4.1: Fixed to keep task requests when no worker is available
    10.19.2018 1.4.4: Worker can use different queue name (=pool) when user task json has 'pool' key
    10.22.2018 1.4.5: Updated workers table; workerId2 for workers table; life_left;
    10.23.2018 1.4.6: Bug fix about workerId2
    10.24.2018 1.4.7: Updated to set -1 for dead workers
    10.26.2018 1.4.8: Fixed user termination error code update (-4)

    10.26.2018 1.5.0: Demo version with static workers tested
    10.28.2018 1.5.1: Added ipaddress to worker hb
    11.02.2018 1.5.2: Updated 'lifeleft' in workers table for the last dead worker

    11.06.2018 1.6.0: Tested static workers and sbatch on cori and denovo
               1.6.1: Remove user account name from queue name for EC2
    11.08.2018 1.6.2: Added custom queue name postfix for testing in Config
    11.13.2018 1.6.3: Updated jtm_submit for large node sbatch; Updated send_hb to use CNAME as
                      postfix; Added jtm_check_manager cli; Added jtm-check-worker;

    11.13.2018 1.7.0: Updated sbatch for static worker cloning; removed 'interval' from hb;
                      Added '-r' jgi-task-manager option; Added CLIENT_HB_RECV_INTERVAL = 8
    11.14.2018 1.7.1: Added dynamic worker spawning

    11.15.2018 1.8.0: Dynamic workers; Changed hb header format;
    11.19.2018 1.8.1: Changed basic_consume callback args; Changed process_task_request to use the
                      custom pool name as task queue if -tp is used; Updated routine to get the
                      current live workers and add a feature to get the num workers per pool name;
    11.20.2018 1.8.2: Fixed message unconsumed from inner result queue
    12.13.2018 1.8.2: Changed mysql to gpdb23

    12.17.2018 1.9.0: Updated to use gpdb23; Test the client on denovo

    12.19.2018 2.0.0: Updated to have multiple jtm instances; Updated jtmInterface to determine
                      task request queue and task result queue in run time; Added jtmHostName to
                      workers table;
    12.21.2018 2.0.1: Updated jtm-submit to use 'cl' param;
    12.22.2018 2.0.2: Testing mysql packet exception error;
    01.03.2019 2.0.3: Changed queue name format to hostname.username.customepoolname;
                      Added exception for losing db connection;
    01.04.2019 2.0.4: Changed to auto_delete=True for main task request/result queues;
                      Changed to auto_delete=True for inner task request/result queues;
    01.08.2019 2.0.5: Changed clonecnt update for handle multiple static workers on a same node
    01.09.2019 2.0.6: Fixed nw option to automatically clone n static workers;

    01.14.2019 3.0.0: Tested with Jaws + jtm + jaws account
    01.24.2019 3.0.2: Updated to print error messages to stderr in jtm cli tools;
                      logger.exception("Failed to call %s.
                      Exit code=%s" % (msg.cmd, msg.returncode))
    01.25.2019 3.0.3: Connected jaws custom pool setting to jtm;
                      Read jtm-submit params and start a pool of workers if necessary;

    01.30.2019 3.1.0: Updated worker and manager to handle wid by option params; Updated to support
                      a custom pool creation for Cromwell scatter function so that a pool can be
                      reused for multiple tasks from a scatter;
    02.01.2019 3.1.1: Set 20 for jtm-submit waittime; Removed all sql pool con;

    02.21.2019 3.2.1: Changed the way to count active workers; ==> not working, reverted 0226
    02.25.2019 3.2.2: Changed to exclusive=False for worker hb queue;
    02.26.2019 3.2.3: Added default queues, small, medium, large, and xlarge;
                      Set default queue as "small";
    02.27.2019 3.2.4: Set default queue as "small"
    02.28.2019 3.2.5: Changed client hb -> worker(s) exchange type from fanout to topic;
    03.06.2019 3.2.6: Changed worker cnt routine to get the # activate workers per pool name;
    03.07.2019 3.2.7: Use unique worker id for processing each individual hb from workers;
                      Changed resource log target dir;
    03.21.2019 3.2.10: Added Cori KNL support;
    03.27.2019 3.2.11: Fixed jtm-submit by checking response value is None or not;

    03.28.2019 3.3.0: Fixed scatter support (Something wrong with exchanges);
                      Redid JtmInterface message receiving;
    04.02.2019 3.3.1: Fixed dynamic worker multi-sbatch; For new worker request, set lifeleft=-2;

    04.02.2019 3.4.0: Improving jtm-kill;
    04.03.2019 3.4.1: Improved process_task_kill() for updated runs table for the case of task
                      cancellation;
    04.03.2019 3.4.3: Changed to single poison queue per clust r;

    04.04.2019 4.0.0: Added a thread for checking termination requested task; Added kill exchnage
                      and queue;
    04.08.2019 4.0.2: Changed result recv interval 6->2secs
    04.15.2019 4.0.3: RESULT_RECEIVE_INTERVAL = 0.5, WORKER_INFO_UPDATE_WAIT = 1;
                      Fixed JtmInterface for recv task id;

    04.15.2019 4.1.0: JtmInterface max wait => 50 rounds;
                      Created separate queue per each jtm interface command;
    04.16.2019 4.2.0: Set each jtm-submit creates a temp queue for recv task id;
    04.17.2019 4.2.2: Fixed getting alive #worker routine -> Set slurm id once sbatched and check
                      slurmid != 0 when count the alive or sbatched workers with given pool name;
                      Restrict auto cloning only if workertype ==1 and slurm jobid > 1;
    04.22.2019 4.2.4: Updated select_count_workers_by_poolname_enddate and
                      select_count_workers_by_jtmhostname to count sbatched worekrs;
                      Added nwpn (num workers per node) to jtm-submit and cromwell conf;
                      Updated process_task_request to double check the number of workers needed;

    04.24.2019 5.0.0: Added nwpn; chance all counting alive workers routines;
                      Added the feature of "shared=0" for jtm-worker pool. If shared=0, the pool won't be shared
                      among workflows which use the same pool name;
    04.25.2019 5.0.2: Updated to detect slurm failure;
    04.29.2019 5.0.3: Added cromwell job id in cromwell conf;
                      Added exclusive;
    04.30.2019 5.0.4: Tested exclusive;

    05.03.2019 5.1.0: Change to insert all workers by nwpn -> Set nWorkersPerNode=1 for dynamic worker;
                      Appended serial number to uniq worker id for dynamic worker -> Change the
                      workerid len from 22 to 23;
                      Change hb interval send->2sec, recv->8sec;
                      Tested jtm-worker -wt static -t 00:05:00 -cl cori -nw 4;
                      Even with nwpn=4, cloning is done by node based (b/c of clone count checking);
                      Resource log subdir -> padded string from task id to store resource log
                      ex) tid=228 --> 00/00/02;

    05.10.2019 5.2.0: Created single db connection for worker's hb recv;
    05.15.2019 5.2.1: Improved jtm-kill by sending kill msg to workers only from task_kill_thread();
    05.16.2019 5.2.2: Updated to record resource log file full path in
    05.22.2019 5.2.3: Added jtm-resource-log;
    05.24.2019 5.2.4: Updated jtminterface wait method from time_limit to sleep;
                      Added taskqueue declaration and biding in jtminterface so that jtm-submit
                      requests are maintained for the case where the manager is not available;
                      Keep jtm-submit waittime = 600s;
                      Bug fix: short task ("ls") status update bug fix in
                      update_runs_tid_startdate_by_tid;

    05.28.2019 5.3.0: Cronjob started on cori20;
    05.30.2019 5.3.1: Tested with py3;
                      pika upgraded to 1.0.1;
                      RmqConnectionHB --> remove heartbeat_interval;
                      no_ack --> auto_ack;
                      basic_consume param changed;
                      cPickle is not supported in py3; Need to upgrade it for py3 in JAWS conda env;
                      jtm needs pip install mysqlclient ==> not working ==> downgrade openssl ==>
                      conda install openssl=1.0.2r (will lower python 3.7.1 to 3.7.0);
                      jaws conda needs "conda install -c conda-forge shortuuid,
                      conda install -c conda-forge pika,
                      conda install -c anaconda numpy";

    06.11.2019 5.4.0: Replace time.sleep() in recv_hb_from_worker()
                      with conn.process_data_events(time_limit=interval) to fix lost connection;
                      Replace time.sleep() in the all callbacks with ch._connection.sleep();
    06.12.2019 5.4.1: Still lost connection in recv_hb_from_worker() -> set heartbeat=600,
                      blocked_connection_timeout=600 in RabbitmqConnecion param;
                      ** single db conn open/close in recv_hb_from_workers_thread;
    06.13.2019 5.4.2: Multiprocessing -> threading;
                      recv_result -> multithreaded (default: 10)
    06.25.2019 5.4.4: Revered back to multiprocessing (ref. https://stackoverflow.com/questions/3044580/multiprocessing-vs-threading-python)
    08.13.2019 5.4.5: Removed threads for on_result(); pika upgraded to v1.1.0;
    09.19.2019 5.4.6: Found missing db connection; Open and close db multiple times in recv_hb_from_workers_thread;
                      Delete tag 5.5.0 (= thread test)
                      Released as a stable production version;

"""

from Config import *
from Common import *
from Sql.SqlStmt import *
from Lib.RabbitmqConnection import *
from Lib.DbUtils import *

# import collections
from math import ceil
from Utils.Run import pad_string_path

print("JGI Task Manager, version: %s" % (VERSION))  # VERSION <- Config.py

# --------------------------------------------------------------------------------------------------
# Globals
# --------------------------------------------------------------------------------------------------
g_nTotalWorkers = multiprocessing.Value('i', 0)


# --------------------------------------------------------------------------------------------------
def recv_hb_from_workers_thread(hbQueueName, parentPid, hbeatProc, logDestDir, resourceReport):
    """
    Receive hearbeats from all the workers running
    :param hbQueueName:
    :param parentPid:
    :param hbeatProc:
    :param logDestDir:
    :param resourceReport:
    :return:
    """

    # TODO: change to event driven by ch.basic_consume and ch.start_consuming
    # Remote broker (mq.nersc.gov)
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()

    exchName = JTM_WORKER_HB_EXCH

    ch.exchange_declare(exchange=exchName,
                        exchange_type="direct",
                        durable=False,
                        auto_delete=False)

    # This queue can be declared from workers first.
    try:
        ch.queue_declare(queue=hbQueueName,
                         durable=False,
                         exclusive=False,
                         auto_delete=True)
    except Exception as detail:
        logger.exception("Exception: The queue, %s is in use.", hbQueueName)
        logger.exception("Detail: %s", str(detail))
        sys.exit(1)

    # Queue binding for the queue to receive hb from workers
    ch.queue_bind(exchange=exchName,
                  queue=hbQueueName,
                  routing_key=hbQueueName)

    bCleared = False
    bFound = False
    numWorkerCheck = 0
    interval = CLIENT_HB_RECV_INTERVAL

    # TODO: need to open/close in place to evade from locking. Any other solutions?
    # Note: tried connection pool but not working as expected
    # db = DbSqlLite(dbpath="jtm.db")
    # db = DbSqlMy(db=MYSQL_DB)

    while 1:
        wids = {}
        try:
            if int(PIKA_VER[0]) < 1:  # v0.13.1
                methodFrame, headerFrame, body = ch.basic_get(queue=hbQueueName, no_ack=True)
            else:  # v1.0.1 or higher
                methodFrame, headerFrame, body = ch.basic_get(queue=hbQueueName, auto_ack=True)

            # Get all the hb messages from workers
            # if methodFrame:
            #    for i in range(int(methodFrame.message_count)):
            #        methodFrame, headerFrame, body = ch.basic_get(queue=hbQueueName, auto_ack=True)
        except Exception as detail:
            logger.exception("Exception: Failed to get a message from %s.", hbQueueName)
            logger.exception("Detail: %s", str(detail))
            ch.close()
            conn.close()
            sys.exit(1)

        if body and not bCleared:
            # To clear up any heartbeat messages left in the heartbeat queue.
            for i in range(int(methodFrame.message_count)):

                if int(PIKA_VER[0]) < 1:  # v0.13.1
                    _, _, _ = ch.basic_get(queue=hbQueueName, no_ack=True)
                else:  # v1.0.1 or higher
                    _, _, _ = ch.basic_get(queue=hbQueueName, auto_ack=True)

            bCleared = True

        elif body and bCleared:
            # To deal with the heartbeats from the worker(s). The workers send
            # each hostname and pid. Using "dict", get the number of unique
            # pids of the running workers.
            msgUnzipped = json.loads(zloads(body))
            # type conversion and sort by key
            msgUnzipped = {int(k): v for k, v in msgUnzipped.items()}

            # worker id is used to collect unique rootPid (= num of workers)
            wid = msgUnzipped[HB_MSG["worker_id"]]
            wids[wid] = msgUnzipped

            # NOTE: Workers send it"s hb interval to the client in the msg packet.
            # Set the checking heartbeat as (jtm-worker"s sending heartbeat
            # interval value * intervalIncRate)

            # Worker's hb should be more than one so consume all the hb's
            for i in range(int(methodFrame.message_count)):
                if int(PIKA_VER[0]) < 1:  # v0.13.1
                    methodFrame, headerFrame, body = ch.basic_get(queue=hbQueueName, no_ack=True)
                else:  # v1.0.1 or higher
                    methodFrame, headerFrame, body = ch.basic_get(queue=hbQueueName, auto_ack=True)

                msgUnzipped = json.loads(zloads(body))
                msgUnzipped = {int(k): v for k, v in msgUnzipped.items()}
                wid = msgUnzipped[HB_MSG["worker_id"]]
                wids[wid] = msgUnzipped

            for k, v in wids.iteritems():
            # for k, v in wids.items():  # py3
                taskId = v[HB_MSG["task_id"]]
                rootPid = v[HB_MSG["root_pid"]]
                childPid = v[HB_MSG["child_pid"]]
                workerId = v[HB_MSG["worker_id"]]
                slurmJobId = v[HB_MSG["slurm_jobid"]]
                workerType = v[HB_MSG["worker_type"]]
                endDate = v[HB_MSG["end_date"]]
                lifeLeft = v[HB_MSG["life_left"]]
                memPerNode = v[HB_MSG["mem_per_node"]]
                memPerCore = v[HB_MSG["mem_per_core"]]
                nCores = v[HB_MSG["num_cores"]]
                jobTime = v[HB_MSG["job_time"]]
                cloneTime = v[HB_MSG["clone_time_rate"]]
                hostName = v[HB_MSG["host_name"]]
                jtmHostName = v[HB_MSG["jtm_host_name"]]
                ipAddr = v[HB_MSG["ip_address"]]
                poolName = v[HB_MSG["pool_name"]]
                numWorkersOnThisNode = v[HB_MSG["nwpn"]]  # TODO: now with nwpn, numWorkersOnThisNode = 1 always

                if resourceReport:
                    logger.resource(v)

                if poolName:
                    qName = JTM_INNER_REQUEST_Q + "." + poolName
                else:
                    qName = ""

                # Fixme: bytearray index out of range EXCEPTION with gpdb23
                success = False
                # TODO: still need to do "while" for checking db connection?
                while success is not True:
                    success = True
                    try:
                        db = DbSqlMy(db=MYSQL_DB)
                        # This increases AUTO_INCREMENT field even there is no insertion
                        # db.execute("""INSERT IGNORE INTO workers (workerId, slurmJobId, workerType)
                        #     VALUES ("%(wid)s", %(slurmjid)d, %(wtype)d)
                        # """ % dict(wid=workerId,
                        #            slurmjid=slurmJobId,
                        #            wtype=workerType))
                        bExists = db.selectScalar(JTM_SQL["select_exists_workers_by_workerid"]
                                                  % dict(wid=workerId))

                        if not bExists:
                            db.execute(JTM_SQL["insert_workers_workerid_slurmjobid"]
                                       % dict(wid=workerId,
                                              slurmjid=slurmJobId,
                                              wtype=workerType,
                                              nwpn=numWorkersOnThisNode))
                        else:
                            # TODO: stress test needed! Test if nworkers > 1000
                            # if workerType == WORKER_TYPE["manual"]
                            #     db.execute(JTM_SQL["update_workers_enddate_lifeleft_by_workerid"]
                            #                % dict(wid=workerId, now=endDate))
                            # elif workerType == WORKER_TYPE["static"]:
                            #     pass
                            # elif workerType == WORKER_TYPE["dynamic"]:
                            #     pass
                            # TODO: can only update endDate and lifeLeft after 1st insert
                            db.execute(JTM_SQL["update_workers_enddate_lifeleft_by_workerid"]
                                       % dict(wid=workerId,
                                              now=endDate,
                                              left=lifeLeft,
                                              mpn=memPerNode if workerType != 0 else 0,
                                              mpc=memPerCore if workerType != 0 else 0,
                                              ncores=nCores if workerType != 0 else 0,
                                              jtime=jobTime,
                                              clonerate=cloneTime if workerType != 0 else 0,
                                              hname=hostName,
                                              jtmhname=jtmHostName,
                                              ipaddr=ipAddr,
                                              poolname=qName,
                                              slurmjid=slurmJobId))

                        db.commit()
                        db.close()
                    except Exception as e:
                        logger.critical(e)
                        logger.critical("Failed to update workers table for enddate and lifeleft.")
                        logger.debug("Retry to update workers table for enddate and lifeleft.")
                        # raise
                        success = False
                        # time.sleep(1)
                        conn.process_data_events(time_limit=1.0)

                if taskId > 0 and rootPid != childPid:  # if there is a user process started
                    dateStr = datetime.datetime.now().strftime("%Y-%m-%d")
                    if logDestDir:
                        logDir = os.path.join(logDestDir, "resource")
                    else:
                        logDir = "%s/resource" % (os.getcwd())

                    paddedDir = pad_string_path(taskId, depth=3)  # 226 --> 00/00/02
                    makeDir(os.path.join(logDir, paddedDir))
                    resourceLogFileName = "%s/%s/jtm_resource_%d_%s.log" % (logDir, paddedDir, taskId, dateStr)
                    # logger.debug("resource log file: %s" % (resourceLogFileName))

                    # new header ################
                    # "child_pid": 1,
                    # "clone_time_rate": 2,
                    # "cpu_load": 3,
                    # "end_date": 4,
                    # "host_name": 5,
                    # "ip_address": 6,
                    # "job_time": 7,
                    # "life_left": 8,
                    # "mem_per_core": 9,
                    # "mem_per_node": 10,
                    # "num_cores": 11,
                    # "num_tasks": 12,
                    # "num_workers_on_node": 13,
                    # "perc_mem_used": 14,
                    # "pool_name": 15,
                    # "ret_msg": 16,
                    # "rmem_usage": 17,
                    # "root_pid": 18,
                    # "run_time": 19,
                    # "slurm_jobid": 20,
                    # "task_id": 21,
                    # "vmem_usage": 22,
                    # "worker_id": 23,
                    # "worker_type": 24
                    # "jtm_host_name": 25
                    # "nwpn": 26
                    # - parent_pid: jtm-worker's pid
                    # - user_command_pid: sh process pid for running user command
                    # - task_id: if registered to JTM, valid task_id, if not, 0

                    with open(resourceLogFileName, 'a') as rf:
                        rf.write(",".join([str(i) for i in v.values()]))
                        rf.write('\n')

                    # Update runs table with resource usage for a task
                    # rsc = ",".join([str(i) for i in collections.OrderedDict(sorted(v.items())).values()])
                    # Update runs table
                    # db.execute(JTM_SQL["update_runs_resources_by_tid"]
                    #            % dict(rsc=rsc, tid=taskId))

                    try:
                        # Update tasks table with "running" status == 2 if status is still 0 or 1
                        db = DbSqlMy(db=MYSQL_DB)
                        # db.execute(JTM_SQL["update_runs_status_workerid_by_tid"]
                        #            % dict(sid=TASK_STATUS["running"],
                        #                   tid=taskId,
                        #                   wid=workerId,
                        #                   cpid=childPid))
                        # db.execute(JTM_SQL["insert_ignore_workers_by_wid_slurmjid"]
                        #            % dict(wid=workerId,
                        #                   slurmjid=slurmJobId))

                        # Update runs table for a task
                        # TODO: really need to store full path to the log?
                        # logger.debug(db.selectAll("select * from runs where taskId=%d" % taskId))
                        db.execute(JTM_SQL["update_runs_status_by_taskid"]
                                   % dict(sid=TASK_STATUS["running"],
                                          tid=taskId,
                                          wid=workerId,
                                          cpid=childPid,
                                          resourcelog=resourceLogFileName))
                        db.commit()
                        db.close()
                    except Exception as e:
                        logger.critical(e)
                        logger.critical("Failed to update runs table for status.")
                        raise

                ####################################################################################
                # Check if cloning is needed
                ####################################################################################
                def hms_to_m(s):
                    t = 0
                    for u in s.split(":"):
                        t = 60 * t + int(u)
                    return int(t / 60)

                jtMin = 0
                if jobTime:
                    jtMin = hms_to_m(jobTime)

                rate = float(lifeLeft) / jtMin if lifeLeft > 0 else 0.0

                if workerType == 1:  # only for static workers
                    logger.debug("worker type {}, wid {}, jobtime {}, jobtime in minute {}, Life left in minute {}, computedRate {}, cloneTimeRate {}".format(
                                  workerType, workerId, jobTime, jtMin, lifeLeft, rate, cloneTime))

                # Static workers cloning is determined by the lifeleft
                # if jobTime and workerType > 0 and slurmJobId > 0 and rate <= float(cloneTime):
                if jobTime and workerType == 1 and slurmJobId > 1 and rate <= float(cloneTime):
                    db = DbSqlMy(db=MYSQL_DB)
                    cloneCnt = db.selectScalar(JTM_SQL["select_clonecnt_workers_by_workerid"]
                                               % dict(wid=workerId))
                    db.close()

                    logger.debug("Clone count = %s" % cloneCnt)

                    if int(cloneCnt) == 0:
                        try:
                            db = DbSqlMy(db=MYSQL_DB)
                            # db.execute(JTM_SQL["update_workers_clonecnt_by_workerid"]
                            #            % dict(wid=workerId))

                            # To handle multiple static workers on a same node (-nw option)
                            # increase clonecnt by job id 01092019
                            db.execute(JTM_SQL["update_workers_clonecnt_by_slurmjobid"]
                                       % dict(slurmjobid=slurmJobId))
                            db.commit()
                            db.close()
                        except Exception as e:
                            logger.critical(e)
                            logger.critical("Failed to update workers table for clonecnt.")
                            raise

                        logger.info("Send cloning signal to the worker, {}".format(workerId))

                        # TODO: if multiple workers are in a node, need to send this reproduce command
                        #  to only one of those workers
                        #  only static worker spawn itself automatically
                        if workerType == WORKER_TYPE["static"]:
                            send_reproduce_or_die(workerId, 0, 1)

                # Dynamic workers cloning and terminating are determined by
                # total # tasks in the queue - number_workers_alive + workers_requested
                #
                # 11.14.2018
                # no automatic cloning of dynamic workers for now
                # user specified custom or base pool will be created and the pool will be disappeared after jobtime
                #
                # if workerType == WORKER_TYPE["dynamic"]:
                #     # TODO: need num clones option, default: 2
                #     # TODO: need to take #tasks in the queue into consideration
                #     send_reproduce_or_die(workerId, 0, 2)
                ########################################################################################################

            #
            # Set unresponsive workers as dead
            #
            # Collect worker_ids from hb
            activeWorkers = []
            for k, v in wids.iteritems():
                activeWorkers.append(v[HB_MSG["worker_id"]])
            # logger.debug("# nique worker IDs in HB = %d" % len(activeWorkers))

            # Collect worker_id which are still set as alive from workers table
            success = False
            while success is not True:
                success = True
                try:
                    db = DbSqlMy(db=MYSQL_DB)
                    # liveWorkers = db.selectAll(JTM_SQL["select_workerid_workers_by_lifeleft"]
                    liveWorkers = db.selectAll(JTM_SQL["select_workerid_workers_by_lifeleft_jtmhostname"]
                                               % dict(jtmhostname=jtmHostName))
                    liveWorkers = [i[0] for i in liveWorkers]
                    # logger.debug("# of unique worker IDs in table = %d" % len(liveWorkers))

                    # Check & update workers table
                    # set lifeLeft as -1 for dead workers
                    for w in liveWorkers:
                        if w not in activeWorkers:
                            db.execute(JTM_SQL["update_workers_lifeleft_by_workerid"]
                                       % dict(wid=w,
                                              jtmhostname=jtmHostName))
                    db.commit()
                    db.close()
                except Exception as e:
                    logger.critical(e)
                    logger.critical("Failed to update workers table for lifeleft.")
                    logger.debug("Retry to update workers table for lifeleft.")
                    # raise
                    success = False
                    # time.sleep(1)
                    conn.process_data_events(time_limit=1.0)

            # g_nTotalWorkers.value = len(activeWorkers)
            db = DbSqlMy(db=MYSQL_DB)
            # liveTotalWorkers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_lifeleft"])
            liveTotalWorkers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_lifeleftt_jtmhostname"]
                                               % dict(jtmhostname=jtmHostName))

            db.close()
            g_nTotalWorkers.value = int(liveTotalWorkers) if liveTotalWorkers else 0
            # logger.debug("Total number of workers in table (alive + requested): %d", g_nTotalWorkers.value)
            logger.debug("# workers: in hb=%d, in table=%d, alive+requested=g_nTotalWorkers=%d" % (len(activeWorkers), len(liveWorkers), g_nTotalWorkers.value))

            if g_nTotalWorkers.value > 0:
                bFound = True
                numWorkerCheck = 0  # reinitialize

        elif not body and bCleared and bFound:
            # If we are here, unfortunately, we lost all the workers that we've been using.
            numWorkerCheck += 1
            g_nTotalWorkers.value = 0
            logger.info("Waiting for worker(s)...")

            # NOTE: 2013.09.05 To prevent from failing to detect workers
            # Increase interval. default=1.2
            # TODO: still need this? ==> 11.13.2018 removed
            # intervalIncRate = intervalIncRate * CLIENT_HB_RECEIVE_INT_INC_RATE

            if WORKER_HB_CHECK_MAX_COUNT != 0 and numWorkerCheck > WORKER_HB_CHECK_MAX_COUNT:  # hit the max checking limit
                # Close connection and kill parent and itself
                ch.close()
                conn.close()
                hbeatProc.terminate()
                os.kill(int(parentPid), signal.SIGTERM)
                # sys.exit(1)
                os._exit(1)

            # If there no workers alive after 5 checks, set lifeLeft to -1 for all
            if numWorkerCheck == 3:
                try:
                    db = DbSqlMy(db=MYSQL_DB)
                    db.execute(JTM_SQL["update_workers_lifeleft_for_last"])
                    db.commit()
                    db.close()
                except Exception as e:
                    logger.critical(e)
                    logger.critical("Failed to update workers table for lifeleft.")
                    raise

        # TODO: Need to start with shorter interval like (0.5-1 sec) for about 30 sec
        # from the beginning and then use the user specified interval
        # ==> 11.13.2018 just wait CLIENT_HB_RECV_INTERVAL
        #time.sleep(interval)

        # To fix connection lost
        # Ref) https://github.com/pika/pika/issues/1224
        try:
            conn.process_data_events(time_limit=float(interval))
        except Exception as e:
            logger.critical(e)
            logger.critical("RMQ connection lost.")
            # sys.exit(1)
            os._exit(1)


    # unreachable
    # db.close()
    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def send_hb_to_workers_thread():
    """
    Broadcast heartbeat to all workers
    Ref) http://www.rabbitmq.com/tutorials/tutorial-three-python.html
    """
    # Remote broker (mq.nersc.gov)
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()
    exchName = JTM_CLIENT_HB_EXCH

    ch.exchange_declare(exchange=exchName,
                        # exchange_type="fanout",
                        exchange_type="topic",
                        durable=False,
                        auto_delete=False)

    msgContainer = {}
    msgContainer["task_type"] = TASK_TYPE["hb"]
    msgZipped = zdumps(json.dumps(msgContainer))

    try:
        while 1:
            ch.basic_publish(exchange=exchName,
                             # routing_key='',
                             routing_key="*." + CNAME,  # all workers with CNAME can hear it
                             body=msgZipped)
            # time.sleep(CLIENT_HB_SEND_INTERVAL)
            conn.process_data_events(time_limit=CLIENT_HB_SEND_INTERVAL)
    except Exception as e:
        logger.critical("Something wrong in send_hb_to_workers_thread(): %s", e)
        raise

    # unreachable
    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def recv_result_from_workers_thread():
    """
    Receive hearbeats from all the workers running
    """
    # Remote broker (mq.nersc.gov)
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()

    exchName = JTM_INNER_MAIN_EXCH
    innerResQName = JTM_INNER_RESULT_Q

    # Default; exchName = jgi_jtm_inner_main_exchange
    ch.exchange_declare(exchange=exchName,
                        exchange_type="direct",
                        passive=False,
                        durable=True,
                        auto_delete=False)
    ch.exchange_declare(exchange=JGI_JTM_MAIN_EXCH,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)

    # This queue can be declared from workers first
    try:
        ch.queue_declare(queue=innerResQName,
                         durable=True,
                         exclusive=False,
                         auto_delete=True)
                         # auto_delete=False)

    except Exception as detail:
        logger.exception("Exception: The queue, %s is in use.", JTM_TASK_RESULT_Q)
        logger.exception("Detail: %s", str(detail))
        sys.exit(1)

    # Queue binding for the queue to receive hb from workers
    ch.queue_bind(exchange=exchName,
                  queue=innerResQName,
                  routing_key=innerResQName)

    # Todo: change to a threaded version
    ch.basic_qos(prefetch_count=NUM_RESULT_RECV_THREADS)

    # NOTE: the below methods might cause error in consuming messages which are already queued
    # ch.basic_consume(lambda ch, method, properties, body: recv_result_on_result(ch, method, properties, body, dbpool),
    #                  queue=innerResQName,
    #                  auto_ack=False)
    # recv_result_on_result_callback = functools.partial(recv_result_on_result, args=(dbpool))
    # ch.basic_consume(recv_result_on_result_callback, innerResQName, auto_ack=False)

    # OLD
    if int(PIKA_VER[0]) < 1:  # v0.13.1
        ch.basic_consume(recv_result_on_result, queue=innerResQName, no_ack=False)
    else:  # v1.0.1 or higher
        ch.basic_consume(queue=innerResQName, on_message_callback=recv_result_on_result)

    # NEW
    # threads = []
    # on_result_callback = functools.partial(on_result, args=(conn, threads))
    # ch.basic_consume(innerResQName, on_result_callback)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        ch.stop_consuming()

    # Wait for all to complete
    # for thread in threads:
    #     thread.join()

    # Fixme: even there are a result in _jtm_inner_result_queue.* queue, it is not consumed.
    # while 1:
    #     methodFrame, headerFrame, body = ch.basic_get(queue=innerResQName, auto_ack=False)
    #     if body:
    #         msgUnzipped = json.loads(zloads(body))
    #         dFlag = int(msgUnzipped["done_flag"])
    #         taskId = int(msgUnzipped["task_id"])
    #         retMsg = msgUnzipped["ret_msg"]
    #         workerId = msgUnzipped["worker_id"]
    #         hostName = msgUnzipped["host_name"]
    #
    #         logger.debug("Result received: {}".format(msgUnzipped))

    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def nack_message(ch, deliveryTag):
    """Note that `ch` must be the same pika channel instance via which
    the message being ACKed was retrieved (AMQP protocol constraint).
    """
    if ch.is_open:
        # ch.basic_ack(delivery_tag)
        ch.basic_reject(delivery_tag=deliveryTag, requeue=False)
    else:
        # Channel is already closed, so we can't ACK this message;
        # log and/or do something that makes sense for your app in this case.
        pass


# -------------------------------------------------------------------------------
def ack_message(ch, deliveryTag):
    """Note that `ch` must be the same pika channel instance via which
    the message being ACKed was retrieved (AMQP protocol constraint).
    """
    if ch.is_open:
        ch.basic_ack(deliveryTag)
    else:
        # Channel is already closed, so we can't ACK this message;
        # log and/or do something that makes sense for your app in this case.
        pass


# -------------------------------------------------------------------------------
# def recv_result_on_result2(conn, ch, delivery_tag, body):
#     # thread_id = threading.get_ident()
#     # LOGGER.info('Thread id: %s Delivery tag: %s Message body: %s', thread_id,
#     #             delivery_tag, body)
#
#     msgUnzipped = json.loads(zloads(body))
#     dFlag = int(msgUnzipped["done_flag"])
#     taskId = int(msgUnzipped["task_id"])
#     retMsg = msgUnzipped["ret_msg"]
#     workerId = msgUnzipped["worker_id"]
#     hostName = msgUnzipped["host_name"]
#
#     logger.info("Result received: {}".format(msgUnzipped))
#
#     if retMsg != "hb":
#         db = DbSqlMy(db=MYSQL_DB)
#         logger.debug("Update tasks {} with {}".format(taskId, dFlag))
#         try:
#             # db = DbSqlMy(db=MYSQL_DB)
#             db.execute(JTM_SQL["update_tasks_doneflag_by_taskid"]
#                        % dict(tid=taskId,
#                               dflag=dFlag))
#             db.commit()
#             # db.close()
#         except Exception as e:
#             logger.critical(e)
#             logger.critical("Failed to update tasks table for doneflag.")
#             raise
#
#         if dFlag > 0:
#             tStatus = TASK_STATUS["success"]  # 4
#         else:
#             if dFlag == -1:
#                 tStatus = TASK_STATUS["outputerror"]
#             elif dFlag == -2:
#                 tStatus = TASK_STATUS["failed"]
#             elif dFlag == -3:
#                 tStatus = TASK_STATUS["outofresource"]
#             elif dFlag == -4:
#                 tStatus = TASK_STATUS["terminated"]
#             else:
#                 logger.critical("Unknown return code {}".format(dFlag))
#                 sys.exit(1)
#
#         logger.debug("Update runs for taskid {} with {}".format(taskId, tStatus))
#
#         workerId2 = 0
#         # db = DbSqlMy(db=MYSQL_DB)
#         rows = db.selectAll(JTM_SQL["select_workerid2_workers_by_wid"]
#                             % dict(wid=workerId))
#         # db.close()
#         try:
#             workerId2 = int(rows[0][0])
#         except:
#             workerId2 = 0
#
#         # Just in case
#         while workerId2 == 0:
#             # db = DbSqlMy(db=MYSQL_DB)
#             rows = db.selectAll(JTM_SQL["select_workerid2_workers_by_wid"]
#                                 % dict(wid=workerId))
#             # db.close()
#             try:
#                 workerId2 = int(rows[0][0])
#             except:
#                 workerId2 = 0
#                 logger.debug(
#                     "select_workerid2_workers_by_wid sleep for %f" % WORKER_INFO_UPDATE_WAIT)
#             # time.sleep(WORKER_INFO_UPDATE_WAIT)
#             ch._connection.sleep(WORKER_INFO_UPDATE_WAIT)
#
#         success = False
#         while success is not True:
#             success = True
#             try:
#                 # db = DbSqlMy(db=MYSQL_DB)
#                 db.execute(JTM_SQL["update_runs_status_workerid2_by_taskid_2"]
#                            % dict(sid=tStatus,
#                                   wid2=workerId2,
#                                   now=time.strftime("%Y-%m-%d %H:%M:%S"),
#                                   tid=taskId))
#                 db.commit()
#                 # logger.debug(db.selectAll("select * from runs where taskId=%d" % taskId))
#                 # db.close()
#
#             except Exception as e:
#                 logger.critical(e)
#                 logger.critical("Failed to update runs table for status and workerid2.")
#                 logger.debug("Retry to update runs table for status and workerid2.")
#                 # raise
#                 success = False
#                 logger.debug(
#                     "update_runs_status_workerid2_by_taskid_2 sleep for %d" % RUNS_INFO_UPDATE_WAIT)
#                 # time.sleep(RUNS_INFO_UPDATE_WAIT)
#                 ch._connection.sleep(RUNS_INFO_UPDATE_WAIT)
#
#         # Print report
#         if dFlag == DONE_FLAGS["success"]:  # 1
#             logger.info("Task %s --> Success on worker/host, %s/%s",
#                         taskId, workerId, hostName)
#
#         elif dFlag == DONE_FLAGS["success with correct output file(s)"]:  # 2
#             logger.info("Task %s --> Success with valid output(s) on worker/host, %s/%s",
#                         taskId, workerId, hostName)
#
#         elif dFlag == DONE_FLAGS["failed to check output file(s)"]:  # -1
#             logger.info("Task %s --> %s, worker/host: %s/%s",
#                         taskId, retMsg, workerId, hostName)
#
#         elif dFlag == DONE_FLAGS["failed to run user command"]:  # -2
#             logger.info(
#                 "Task %s --> Failed with non-zero exit code. stdout = %s, worker/host: %s/%s",
#                 taskId, retMsg, workerId, hostName)
#
#         elif dFlag == DONE_FLAGS["failed with out-of-mem"]:  # -3
#             pass
#
#         elif dFlag == DONE_FLAGS["failed with user termination"]:  # -4
#             logger.info("Task %s --> Failed by user termination. stdout = %s, worker/host: %s/%s",
#                         taskId, retMsg, workerId, hostName)
#         else:
#             logger.warning("Cannot recognize the return code: %d" % dFlag)
#             logger.info("Task %s --> worker/host: %s/%s", taskId, workerId, hostName)
#             sys.exit(1)
#
#         db.close()
#
#     else:
#         # If it is not a result msg, return it back to the exchange
#         # ch.basic_reject(delivery_tag=delivery_tag, requeue=False)
#         cb2 = functools.partial(nack_message, ch, delivery_tag)
#         conn.add_callback_threadsafe(cb2)
#
#     # After this the result message will be deleted from RabbitMQ
#     # If this worker crashes while running a user command, this task will
#     # be sent to other workers available
#     # ch.basic_ack(delivery_tag=delivery_tag)
#     cb = functools.partial(ack_message, ch, delivery_tag)
#     conn.add_callback_threadsafe(cb)


# -------------------------------------------------------------------------------
# def on_result(ch, method_frame, _header_frame, body, args):
#     (conn, thrds) = args
#     delivery_tag = method_frame.delivery_tag
#     t = threading.Thread(target=recv_result_on_result2, args=(conn, ch, delivery_tag, body))
#     t.start()
#     thrds.append(t)


# -------------------------------------------------------------------------------
def recv_result_on_result(ch, method, props, body):
    """
    recv_result_from_workers_thread's basic_consume callback
    :param ch: channel
    :param method: methodFrame
    :param props: pika connection property
    :param body: message received
    :return:
    """
    if body:
        msgUnzipped = json.loads(zloads(body))
        dFlag = int(msgUnzipped["done_flag"])
        taskId = int(msgUnzipped["task_id"])
        retMsg = msgUnzipped["ret_msg"]
        workerId = msgUnzipped["worker_id"]
        hostName = msgUnzipped["host_name"]

        # logger.debug("Result received: {} {}".format(method, props))
        logger.info("Result received: {}".format(msgUnzipped))

        if retMsg != "hb":
            logger.debug("Update tasks {} with {}".format(taskId, dFlag))
            try:
                db = DbSqlMy(db=MYSQL_DB)
                db.execute(JTM_SQL["update_tasks_doneflag_by_taskid"]
                           % dict(tid=taskId,
                                  dflag=dFlag))
                db.commit()
                db.close()
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed to update tasks table for doneflag.")
                raise

            # Check the task status code
            if dFlag > 0:
                tStatus = TASK_STATUS["success"]  # 4
            else:
                if dFlag == -1:
                    tStatus = TASK_STATUS["outputerror"]
                elif dFlag == -2:
                    tStatus = TASK_STATUS["failed"]
                elif dFlag == -3:
                    tStatus = TASK_STATUS["outofresource"]
                elif dFlag == -4:
                    tStatus = TASK_STATUS["terminated"]
                else:
                    logger.critical("Unknown return code {}".format(dFlag))
                    sys.exit(1)

            # This seems like resolving runs table lock issue
            # issue: status is not changed to 4 after the update
            # TODO: need to improve
            #############################################
            # time.sleep(RESULT_RECEIVE_INTERVAL)
            ch._connection.sleep(RESULT_RECEIVE_INTERVAL)
            #############################################

            logger.debug("Update runs for taskid {} with {}".format(taskId, tStatus))

            # Sometimes workerId2 ==> 0
            # so wait a little bit if that happened
            workerId2 = 0
            db = DbSqlMy(db=MYSQL_DB)
            rows = db.selectAll(JTM_SQL["select_workerid2_workers_by_wid"]
                                % dict(wid=workerId))
            db.close()
            try:
                workerId2 = int(rows[0][0])
            except:
                workerId2 = 0

            # Just in case
            while workerId2 == 0:
                db = DbSqlMy(db=MYSQL_DB)
                rows = db.selectAll(JTM_SQL["select_workerid2_workers_by_wid"]
                                    % dict(wid=workerId))
                db.close()
                try:
                    workerId2 = int(rows[0][0])
                except:
                    workerId2 = 0
                    # logger.debug("select_workerid2_workers_by_wid sleep for %d" % WORKER_INFO_UPDATE_WAIT)

                # time.sleep(WORKER_INFO_UPDATE_WAIT)
                ch._connection.sleep(WORKER_INFO_UPDATE_WAIT)

            # new
            # Fixme: bytearray index out of range EXCEPTION from gpdb23. Seems like network delay to
            #   the mysql server
            # note: 01032019 fixed by while for checking success
            # Fixme: still seeing the bytearray index out of range EXCEPTION -> revert back to non-pool
            success = False
            while success is not True:
                success = True
                try:
                    db = DbSqlMy(db=MYSQL_DB)
                    db.execute(JTM_SQL["update_runs_status_workerid2_by_taskid_2"]
                               % dict(sid=tStatus,
                                      wid2=workerId2,
                                      now=time.strftime("%Y-%m-%d %H:%M:%S"),
                                      tid=taskId))
                    db.commit()
                    # logger.debug(db.selectAll("select * from runs where taskId=%d" % taskId))
                    db.close()

                except Exception as e:
                    logger.critical(e)
                    logger.critical("Failed to update runs table for status and workerid2.")
                    logger.debug("Retry to update runs table for status and workerid2.")
                    # raise
                    success = False
                    logger.debug("update_runs_status_workerid2_by_taskid_2 sleep for %d" % RUNS_INFO_UPDATE_WAIT)
                    # time.sleep(RUNS_INFO_UPDATE_WAIT)
                    ch._connection.sleep(RUNS_INFO_UPDATE_WAIT)

            # Print report
            if dFlag == DONE_FLAGS["success"]:  # 1
                logger.info("Task %s --> Success on worker/host, %s/%s",
                            taskId, workerId, hostName)

            elif dFlag == DONE_FLAGS["success with correct output file(s)"]:  # 2
                logger.info("Task %s --> Success with valid output(s) on worker/host, %s/%s",
                            taskId, workerId, hostName)

            elif dFlag == DONE_FLAGS["failed to check output file(s)"]:  # -1
                logger.info("Task %s --> %s, worker/host: %s/%s",
                            taskId, retMsg, workerId, hostName)

            elif dFlag == DONE_FLAGS["failed to run user command"]:  # -2
                logger.info("Task %s --> Failed with non-zero exit code. stdout = %s, worker/host: %s/%s",
                            taskId, retMsg, workerId, hostName)

            elif dFlag == DONE_FLAGS["failed with out-of-mem"]:  # -3
                pass

            elif dFlag == DONE_FLAGS["failed with user termination"]:  # -4
                logger.info("Task %s --> Failed by user termination. stdout = %s, worker/host: %s/%s",
                            taskId, retMsg, workerId, hostName)
            else:
                logger.warning("Cannot recognize the return code: %d" % dFlag)
                logger.info("Task %s --> worker/host: %s/%s", taskId, workerId, hostName)
                sys.exit(1)

        else:
            # If it is not a result msg, return it back to the exchange
            ch.basic_reject(delivery_tag=method.delivery_tag, requeue=False)

        # After this the result message will be deleted from RabbitMQ
        # If this worker crashes while running a user command, this task will
        # be sent to other workers available
        ch.basic_ack(delivery_tag=method.delivery_tag)

    else:
        logger.critical("No data found from the result message.")
        sys.exit(1)


# -------------------------------------------------------------------------------
def process_task_request(ch, method, props, msg, innerTaskReqQ):
    """
    Get task request from jtm-submit and send it to a worker
    :param ch:
    :param method:
    :param props:
    :param msg: unzipped dict
    :param innerTaskReqQ: inner task queue name (jgi-task-manager --> worker)
    :return:
    """
    # Parse the msg from jtm-submit
    # Example task json
    # {
    #     "command": "ps -aef | grep jtm > /tmp/jtm_ps_grep.out",
    #     "output_dir": "/tmp",
    #     "output_files": "/tmp/jtm_ps_grep.out",
    #     "pool": {"name": "test",
    #              "node": 1,
    #              "nwpn": 1,
    #              "cluster": "cori",
    #              "time": "00:10:00",
    #              "cpu": 1,
    #              "mem": "1GB"}
    # }
    if "command" not in msg:
        logger.critical("Critical: cannot find user command.")
        send_msg_callback(ch, method, props, -5, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)

    if "task_type" not in msg:
        logger.critical("Critical: cannot find task type.")
        send_msg_callback(ch, method, props, -5, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)

    userTask = msg["command"]
    taskType = msg["task_type"]
    outFile = msg["output_files"] if "output_files" in msg else ""  # comma separated list ex) "a.out,b.out,c.out"
    outDir = msg["output_dir"] if "output_dir" in msg else ""
    stdoutFile = msg["stdout"] if "stdout" in msg else ""
    stderrFile = msg["stderr"] if "stderr" in msg else ""
    # cromwellJid = msg["job_id"] if "job_id" in msg else ""

    # Deal with custom pool of workers
    # If user task json has "pool", jtm creates the specified number (=size) of workers
    # and processes the user task using the workers.
    # pool.name: new custom pool name
    # pool.size: number of workers to spawn in the pool
    # pool.cluster: HPC name to use
    # pool.time: worker runtime
    # pool.ncpus: each worker's ncpus requirement (for slurm)
    # pool.mem: each worker' mem requirement (for slurm)
    #
    # ex)
    # $ jtm-submit -cr 'ls' -cl cori -p test2 -t "00:10:00" -c 1 -s 1 -m 5G
    #
    bFailedToRequestWorker = False

    if "pool" in msg and "name" in msg["pool"] and "time" in msg["pool"]:
        # {u'resource': u'cori', u'name': u'test', u'size': 1}
        # ==> {"resource": "cori", "name": "test", "size": 1}
        poolJson = json.loads(json.dumps(msg["pool"]))

        # Set default values defined in Config.py
        # poolName = poolJson["name"] if "name" in poolJson and poolJson["name"] else None
        poolName = poolJson["name"]
        poolCluster = CLUSTER
        # poolTime = poolJson["time"]
        poolNcpus = NCPUS
        poolMem = MEMPERNODE
        poolConstraint = CORI_CONSTRAINT
        poolQos = CORI_QOS
        poolAccount = CORI_ACCNT

        # Note: pool size = numNodes * numWorkersPerNode
        numWorkersPerNode = NWORKERS
        numNodes = NNODES

        # Worker type is restricted to "dynamic" for now.
        # TODO: add the feature to remove the custom pool by user task to support "static"
        #   workers in custom pool
        # workerType = "dynamic"

        # NOTE: if pool_size is set, which means user sets the number of workers needed,
        # n number of jtm dynamic workers will be created, and the pool of n dynamic
        # workers won't be increased over n. The workers will be terminated if there is
        # no tasks requested
        #
        # if not,
        # only one dynamic worker will be created. The worker will be terminated if there is
        # no tasks in the queue for a specified time duration

        # if "size" in poolJson:  # pool size = the number workers in the pool
        #     poolSize = int(poolJson["size"])
        if "cluster" in poolJson:  # cluster/clouds name
            poolCluster = poolJson["cluster"]
        if "time" in poolJson:  # wallclocktime request
            poolTime = poolJson["time"]
        if "cpu" in poolJson:  # number of cores request
            poolNcpus = int(poolJson["cpu"])
        if "mem" in poolJson:  # node memory request
            poolMem = poolJson["mem"]
        if "constraint" in poolJson:  # [haswell | knl | skylake]
            poolConstraint = poolJson["constraint"]
        if "qos" in poolJson:
            poolQos = poolJson["qos"]  # ["genepool_special", "genepool_shared", "jgi_shared", "jgi_exvivo"]
        if "shared" in poolJson:
            # Note: This is not used for now.
            #  If shared=0 is set from jtm-submit and if the pool needs to be terminated forcefully,
            #  jtm-manager should send a poison to the worker(s).
            poolShared = poolJson["shared"]
        if "nwpn" in poolJson:
            numWorkersPerNode = int(poolJson["nwpn"])
        if "node" in poolJson:
            numNodes = int(poolJson["node"])

        assert (len(userTask) <= 1024)
        assert (len(outFile) <= 1024)

        # Create pool
        # Step
        # 1. check the number if workers in the custom pool
        # 2. if 0, create a pool
        #    else send tasks to the pool

        # TODO: maintain the number of workers sbatched -->
        #   nWorkerNeeded = poolSize - nLiveWorkerInPool - nWorkerSbatched
        db = DbSqlMy(db=MYSQL_DB)
        # nLiveWorkerInPool = db.selectScalar(JTM_SQL["select_count_workers_by_poolname_enddate"]
        #                                     % dict(poolname=innerTaskReqQ,
        #                                            hbinterval=WORKER_HB_RECV_INTERVAL * 2))

        # Note: This is also to check the number of sbatched workers in the case of scatter operation in Cromwell
        # If the # of sbatched dynamic worker(s) is less than the requested pool size,
        # then do sbatch for the rest (=nWorkerNeeded)
        #
        # nLiveWorkerInPool = db.selectScalar(JTM_SQL["select_count_workers_by_poolname"]
        #                                     % dict(poolname=innerTaskReqQ,
        #                                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))

        nLiveNodesInPool = db.selectScalar(JTM_SQL["select_count_distinct_jid_workers_by_poolname"]
                                           % dict(poolname=innerTaskReqQ,
                                                  hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("select_count_distinct_jid_workers_by_poolname: %s"
                     % JTM_SQL["select_count_distinct_jid_workers_by_poolname"]
                     % dict(poolname=innerTaskReqQ,
                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        db.close()

        poolSize = numNodes * numWorkersPerNode
        nLiveWorkerInPool = nLiveNodesInPool * numWorkersPerNode
        nWorkerNeeded = int(ceil(float(poolSize - nLiveWorkerInPool) / numWorkersPerNode))

        logger.debug("nWorkerNeeded={} poolSize={} nLiveWorkerInPool={} nLiveNodesInPool={} g_nTotalWorkers={}".format(
                      nWorkerNeeded, poolSize, nLiveWorkerInPool, nLiveNodesInPool, g_nTotalWorkers.value))

        uniqWorkerId = None
        for i in range(0, nWorkerNeeded):
            bFailedToRequestWorker = False
            # NOTE: if poolName is None, JTM will spawn worker(s) in the main pool, not in a custom pool
            # sbatchCmd = """ssh -t {}@{}.nersc.gov "{} && jtm-worker -wt dynamic {} -cl {} -c {} -t {} -m {}" """
            # Now only consider jtm running from nersc nodes and multiple jtm instances per each platform
            #
            # NOTE: User can request only "dynamic" workers from WDL. The "static" workers are managed
            #  by the admin.
            uniqWorkerId = str(shortuuid.uuid())
            sbatchCmd = """{} && jtm-worker -wt dynamic {} -cl {} -c {} -t {} -m {} -wi {} -C {} -nw {}""".format(
                           ENV_ACTIVATION,
                           "-tp %s" % poolName if poolName else "",
                           poolCluster,
                           poolNcpus,
                           poolTime,
                           poolMem,
                           uniqWorkerId,
                           poolConstraint,
                           numWorkersPerNode)

            logger.info("Executing {}".format(sbatchCmd))
            so, se, ec = run_sh_command(sbatchCmd, live=True, log=logger)

            # Print job script for logging
            run_sh_command(sbatchCmd + " --dry-run", live=True, log=logger)

            # Get the slurm job id returned from jtm-worker
            try:
                slurmJobId = int(so.split('\n')[1])
            except:
                logger.critical("Failed to get a valid job ID back from requesting a dynamic worker")
                ec = 1  # make it fail

            if ec == 0:  # if sbatch by jtm-worker done successfully
                logger.debug("Insert into workers table.")
                for nwpn in range(numWorkersPerNode):
                    # Insert into workers for new worker id
                    try:
                        db = DbSqlMy(db=MYSQL_DB)
                        # If a dynamic worker is requested successfully,
                        # insert the info into workers table
                        # loggger.info("Try to update workers table")check_worker

                        # Fixme: After sbatch, another sbatch with the same pool name will be executed
                        #  again.
                        # Solution: Set the lifeleft=-2 and update select_count_workers_by_poolname sql
                        # statement to check only lifeleft!=-1 so that nLiveWorkerInPool can include
                        # already sbatched workers for the pool.
                        #
                        db.execute(JTM_SQL["insert_workers_workerid_workertype_poolname"]
                                   % dict(wid=uniqWorkerId + str(nwpn+1),
                                          wtype=WORKER_TYPE["dynamic"],
                                          poolname=innerTaskReqQ,
                                          jtmhostname=poolCluster,
                                          lifeleft=-2,
                                          slurmjobid=slurmJobId,
                                          # nwpn=numWorkersPerNode,
                                          nwpn=1))
                        db.commit()
                        db.close()
                    except Exception as e:
                        logger.critical(e)
                        logger.critical("Failed to insert workers table for workerid and workertype.")
                        logger.debug("Retry to insert workers table for workerid and workertype.")
                        raise
            else:
                logger.critical("Failed to execute the command, %s" % (sbatchCmd))
                logger.critical("Failed to request workers.")
                send_msg_callback(ch, method, props, TASK_STATUS["invalidtask"], JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, deliveryMode=2)
                bFailedToRequestWorker = True
                break

    if not bFailedToRequestWorker:
        lastTid = -1

        # Fixme: mysql connection.py IndexError: bytearray index out of range
        # note: seems like mysql pool connection error
        success = False
        while success is not True:
            success = True
            try:
                db = DbSqlMy(db=MYSQL_DB)
                # table fields: userCmd, outFiles, doneFlag, retryCnt, taskType
                db.execute(JTM_SQL["insert_tasks_usercmd_outfiles"]
                           % (userTask, outFile, "0", 0, TASK_TYPE[taskType]))
                lastTid = db.selectScalar(JTM_SQL["select_last_insert_id"])
                logger.debug("lastTid = %d" % lastTid)
                db.commit()
                db.close()
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed to insert tasks table for new user task.")
                logger.debug("Retry to insert tasks table for a user task.")
                # lastTid = -1
                # raise  # testing
                success = False
                # time.sleep(1)
                ch._connection.sleep(1.0)

        success = False
        while success is not True:
            success = True
            try:
                db = DbSqlMy(db=MYSQL_DB)
                db.execute(JTM_SQL["insert_runs_tid_sid"]
                           % dict(tid=lastTid,
                                  sid=TASK_STATUS["ready"]))
                db.commit()
                db.close()
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed to insert runs table for a new run.")
                logger.debug("Retry to insert runs table for a new run.")
                # lastTid = -1
                # raise  # testing
                success = False
                # time.sleep(1)
                ch._connection.sleep(1.0)

        # Todo: 01282019 now every cluster has a jtm manager running, so it doesn't need to spawn
        #  dynamic worker by a static worker
        if lastTid != -1:  # if it successfully updates runs table and gets a task id
            # Todo: Check if cancelled or terminated
            db = DbSqlMy(db=MYSQL_DB)
            tStatus = int(db.selectScalar(JTM_SQL["select_status_runs_by_taskid"]
                                          % dict(tid=lastTid)))
            bCancelled = int(db.selectScalar(JTM_SQL["select_cancelled_runs_by_taskid"]
                                             % dict(tid=lastTid)))
            db.close()

            # Note: this is just in case
            #  It is based on the assumption that a task can be cancelled between ready -> queued
            #  status change.
            if bCancelled == 1 or tStatus == TASK_STATUS["terminated"]:
                if tStatus != TASK_STATUS["terminated"]:
                    db = DbSqlMy(db=MYSQL_DB)
                    db.execute(JTM_SQL["update_tasks_doneflag_by_taskid"]
                               % dict(tid=lastTid,
                                      dflag=TASK_STATUS["terminated"]))
                    db.commit()
                    db.close()

            else:
                # Prepare msg to jtm-worker
                msgContainer = {}
                msgContainer["task_id"] = lastTid
                msgContainer["user_cmd"] = userTask
                msgContainer["output_files"] = outFile
                msgContainer["done_flag"] = 0
                msgContainer["task_type"] = TASK_TYPE[taskType]
                msgContainer["output_dir"] = outDir
                msgContainer["stdout"] = stdoutFile
                msgContainer["stderr"] = stderrFile
                # msgContainer["cromwell_jid"] = cromwellJid

                logger.info("Total number of workers (alive + requested): %d", g_nTotalWorkers.value)

                # Create and send request message to workers
                msgZipped = zdumps(json.dumps(msgContainer))
                corrId = str(uuid.uuid4())

                logger.info("Send a task to {}".format(innerTaskReqQ))

                try:
                    ch.basic_publish(exchange=JTM_INNER_MAIN_EXCH,
                                     routing_key=innerTaskReqQ,
                                     properties=pika.BasicProperties(
                                         delivery_mode=2,  # make message persistent
                                         reply_to=JTM_INNER_RESULT_Q,  # set reply queue name
                                         correlation_id=corrId),
                                     body=msgZipped)
                except Exception as detail:
                    logger.exception("Exception: Failed to send a request to %s", innerTaskReqQ)
                    logger.exception("Detail: %s", str(detail))
                    # Todo: set the task status --> failed
                    sys.exit(1)

                # Update status to "queued"
                # TODO: need this update to change the task status to "queued"?
                time.sleep(TASK_STAT_UPDATE_INTERVAL)

                try:
                    db = DbSqlMy(db=MYSQL_DB)
                    db.execute(JTM_SQL["update_runs_tid_startdate_by_tid"]
                               % dict(sid=TASK_STATUS["queued"],
                                      now=time.strftime("%Y-%m-%d %H:%M:%S"),
                                      tid=lastTid))
                    db.commit()
                    db.close()
                except Exception as e:
                    logger.critical(e)
                    logger.critical("Failed to update runs table for status and startdate.")
                    # Todo: properly update runs for the failure
                    # Todo: set the task status --> failed
                    # lastTid = -1
                    raise

        # Send task id to jtm_receive
        logger.debug("Sending task id, %d to jtm-submit via %s" % (lastTid, props.reply_to))
        send_msg_callback(ch, method, props, lastTid, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, deliveryMode=2)


# -------------------------------------------------------------------------------
def process_task_status(ch, method, props, taskId):
    """
    Process task status request
    :param ch: channel
    :param method: methodFrame
    :param props: pika connection property
    :param taskId: for task ID
    :return:
    """
    db = DbSqlMy(db=MYSQL_DB)
    cur = db.execute(JTM_SQL["select_status_runs_by_taskid"] % dict(tid=taskId))
    ret = cur.fetchone()

    if ret:
        taskStatus = ret[0]
    else:
        taskStatus = -1
    db.close()

    # Send tid to jtm_status
    logger.info("Sent reply: %s via %s" % (str(taskStatus), props.reply_to))
    send_msg_callback(ch, method, props, taskStatus, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, deliveryMode=2)


# -------------------------------------------------------------------------------
def process_resource_log(ch, method, props, taskId):
    """
    With given task id, get the resource log file and send to jtm-resource-log

    :param ch:
    :param method:
    :param props:
    :param taskId:
    :return:
    """
    resourceFile = None
    db = DbSqlMy(db=MYSQL_DB)
    cur = db.execute(JTM_SQL["select_resource_runs_by_taskid"] % dict(tid=taskId))
    ret = cur.fetchone()
    db.close()
    if ret:
        resourceFile = ret[0]

    # Send tid to jtm_status
    logger.info("Sent reply: %s via %s" % (resourceFile, props.reply_to))
    send_msg_callback(ch, method, props, resourceFile, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, deliveryMode=2)


# -------------------------------------------------------------------------------
def send_task_kill_request(tid, wid, cpid):
    """
    Send a task kill request
    1. if wid is ready, send wid with tid
    2. else any worker in not ready, so just send tid

    :param tid: task id
    :param wid: worker id
    :param cpid: child process id
    :return:
    """
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()

    exch = JTM_TASK_KILL_EXCH
    qName = JTM_TASK_KILL_Q

    ch.exchange_declare(exchange=exch,
                        exchange_type="fanout",
                        durable=True,
                        auto_delete=False)

    message = {"task_id": tid,
               "worker_id": wid,
               "child_pid": cpid}
    msgZipped = zdumps(json.dumps(message))

    ch.queue_declare(queue=qName,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)

    # Create routing key with "wid"."tid" so that workers can filter the messages
    rKey = str(wid)
    ch.queue_bind(exchange=exch,
                  queue=qName,
                  routing_key=rKey)

    try:
        logger.info("Send send_task_kill_request to worker %s, %r with routing key %s"
                    % (qName, message, rKey))
        assert qName.endswith(CNAME)
        ch.basic_publish(exchange=exch,
                         routing_key=qName,
                         properties=pika.BasicProperties(delivery_mode=2),  # make message persistent
                         body=msgZipped)
    except Exception as e:
        logger.critical("Something wrong in send_task_kill_request(): %s", e)
        raise

    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def send_reproduce_or_die(wid, cpid, nClones):
    """
    Let a worker kill itself or clone itself
     - static worker: nClones = 1
     - dynamic worker: nClones >= 2

    * Note: dynamic workers are no longer cloned. Thus, this function is only for cloning static
            workers.

    :param wid:  worker id
    :param cpid: child process id
    :param nClones: number of clones to spawn (if -1, kill the worker)
    :return:
    """
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()
    ch.exchange_declare(exchange=JTM_WORKER_POISON_EXCH,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)

    qName = JTM_WORKER_POISON_Q
    message = {"worker_id": wid,
               "child_pid": cpid,
               "num_clones": nClones}
    msgZipped = zdumps(json.dumps(message))

    ch.queue_declare(queue=qName,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)
    ch.queue_bind(exchange=JTM_WORKER_POISON_EXCH,
                  queue=qName,
                  routing_key=wid)

    try:
        if nClones == -1:
            logger.info("Send poison to worker %s, %r" % (qName, message))
        else:
            logger.info("Send cloning signal to worker %s, %r" % (qName, message))

        assert qName.endswith(CNAME)
        ch.basic_publish(exchange=JTM_WORKER_POISON_EXCH,
                         routing_key=wid,
                         properties=pika.BasicProperties(delivery_mode=2),  # make message persistent
                         body=msgZipped)
        logger.debug("send_reproduce_or_die published {} with {}".format(message, wid))

    except Exception as e:
        logger.critical("Something wrong in send_poison(): %s", e)
        raise

    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def process_task_kill(ch, method, props, msg):
    """
    Process request to terminate a task
    """
    # Parse the msg from jtm-kill
    # Get workerid
    # Delete tasks row for task_id (runs table will be cascaded)
    # Send kill msg
    taskId = int(msg["task_id"])
    retMsg = None

    db = DbSqlMy(db=MYSQL_DB)
    try:
        tStatus = db.selectScalar(JTM_SQL["select_status_runs_by_taskid"] % dict(tid=taskId))
    except AssertionError:
        tStatus = None
    db.close()


    def update_runs_cancelled(tid):
        db = DbSqlMy(db=MYSQL_DB)
        db.execute(JTM_SQL["update_runs_cancelled_by_tid"]
               % dict(tid=tid,
                      now=time.strftime("%Y-%m-%d %H:%M:%S")))
        db.commit()
        db.close()


    if tStatus:
        # All task status
        # "ready": 0,
        # "queued": 1,
        # "running": 2,
        # "success": 4,
        # "outputerror": -1,
        # "failed": -2,
        # "outofresource": -3,
        # "terminated": -4,
        # "invalidtask": -5
        if tStatus in (TASK_STATUS["ready"], TASK_STATUS["queued"]):
            logger.debug("Task cancellation requested but the task has already been queued. The task will be terminated once it's started.")
            # Update runs table with cancellation requested (cancelled = 1)
            update_runs_cancelled(taskId)
            retMsg = 0
        elif tStatus == TASK_STATUS["running"]:
            logger.debug("Task cancellation requested. The task is being terminated.")
            # # Now worker id and child pid is ready
            # db = DbSqlMy(db=MYSQL_DB)
            # wid = db.selectAs1Col(JTM_SQL["select_workerid_workers_by_tid"] % dict(tid=taskId))
            # childPid = db.selectScalar(JTM_SQL["select_chilpid_runs_by_tid"] % dict(tid=taskId))
            # logger.debug("SQL: select_chilpid_runs_by_tid = %d" % (int(childPid) if childPid else 0))
            # logger.debug("SQL: select_workerid_workers_by_tid %d = %s" % (taskId, wid))
            # db.close()
            #
            # assert childPid > 0
            # assert wid is not None
            #
            # send_task_kill_request(taskId, wid[0], childPid)

            # Update runs table with cancellation requested (cancelled = 1)
            update_runs_cancelled(taskId)
            retMsg = 0
        elif tStatus == TASK_STATUS["success"]:
            logger.debug("Task cancellation requested but the task is in completed status.")
            retMsg = 0
        elif tStatus == TASK_STATUS["terminated"]:
            logger.debug("Task cancellation requested but the task is already in terminated status.")
            retMsg = 0
        elif tStatus == TASK_STATUS["failed"]:
            logger.debug("Task cancellation requested but the task is already in failed status.")
            retMsg = 0
        elif tStatus in (TASK_STATUS["outputerror"], TASK_STATUS["outputerror"], TASK_STATUS["outputerror"]):
            logger.debug("Task cancellation request is ignored.")
            retMsg = 0
        else:
            logger.debug("Failed to cancel a task. Unexpected condition.")
            retMsg = -1
    else:  # task id not found
        retMsg = -5

    # Send status to jtm_kill
    logger.debug("Return kill command result, %d to jtm-kill." % (retMsg))
    send_msg_callback(ch, method, props, retMsg, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, deliveryMode=2)


# -------------------------------------------------------------------------------
def process_check_worker(ch, method, props, msgUnzipped):
    # If custom_pool name is set, get the number of workers in the pool
    # else the total number of live workers will be sent
    logger.debug("jtm-check-worker: %s" % str(msgUnzipped))
    if "task_pool" in msgUnzipped and msgUnzipped["task_pool"]:
        db = DbSqlMy(db=MYSQL_DB)
        # Try to select "poolName" in workers table by hostname + username + poolname
        # Check timediff(now()-endDate) in workers table to filter out dead worker

        nLiveWorkerInPool = db.selectScalar(JTM_SQL["select_count_workers_by_poolname_enddate"]
                                            % dict(poolname=JTM_INNER_REQUEST_Q + '.' + msgUnzipped["task_pool"],
                                                   hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("select_count_workers_by_poolname_enddate: %s"
                     % JTM_SQL["select_count_workers_by_poolname_enddate"]
                     % dict(poolname=JTM_INNER_REQUEST_Q + '.' + msgUnzipped["task_pool"],
                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))

        nTotalWorkers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_poolname_enddate"]
                                        % dict(poolname=JTM_INNER_REQUEST_Q + '.' + msgUnzipped["task_pool"],
                                               hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("select_sum_nwpn_workers_by_poolname_enddate: %s"
                     % JTM_SQL["select_sum_nwpn_workers_by_poolname_enddate"]
                     % dict(poolname=JTM_INNER_REQUEST_Q + '.' + msgUnzipped["task_pool"],
                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        db.close()

        logger.debug("node cnt in the pool: %s" % nLiveWorkerInPool)
        logger.debug("worker cnt in the pool: %s" % nTotalWorkers)

        # send_msg_callback(ch, method, props, nLiveWorkerInPool, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)
        send_msg_callback(ch, method, props, nTotalWorkers if nTotalWorkers else 0, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)

    elif "jtm_host_name" in msgUnzipped and msgUnzipped["jtm_host_name"]:
        db = DbSqlMy(db=MYSQL_DB)
        # Check timediff(now()-endDate) in workers table to filter out dead worker
        nLiveWorkers = db.selectScalar(JTM_SQL["select_count_workers_by_jtm_host_name"]
                                       % dict(jtmhostname=msgUnzipped["jtm_host_name"],
                                              hbinterval=WORKER_HB_RECV_INTERVAL * 3))

        nTotalWorkers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_jtm_host_name"]
                                        % dict(jtmhostname=msgUnzipped["jtm_host_name"],
                                               hbinterval=WORKER_HB_RECV_INTERVAL * 3))

        db.close()
        logger.debug("node cnt in the host: %s" % nLiveWorkers)
        logger.debug("worker cnt in the host: %s" % nTotalWorkers)

        # send_msg_callback(ch, method, props, nLiveWorkers, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)
        send_msg_callback(ch, method, props, nTotalWorkers if nTotalWorkers else 0, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)
    else:
        send_msg_callback(ch, method, props, g_nTotalWorkers.value, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)


# -------------------------------------------------------------------------------
def process_remove_pool(ch, method, props, msgUnzipped):
    taskPool = JTM_INNER_REQUEST_Q + '.' + msgUnzipped["task_pool"]
    db = DbSqlMy(db=MYSQL_DB)
    logger.debug(JTM_SQL["select_jid_workers_by_poolname"] % dict(poolname=taskPool))

    # Get the list of slurm job id to cancel
    zombieJid = db.selectAll(JTM_SQL["select_jid_workers_by_poolname"] % dict(poolname=taskPool))
    for jid in zombieJid:
        # print jid[0]
        scancelCmd = "scancel %s" % (jid[0])
        so, se, ec = run_sh_command(scancelCmd, live=True, log=logger)
        if ec == 0:
            logger.info("Successfully cancel the job, %s" % (jid[0]))
        else:
            logger.debug("%s not found." % (jid[0]))

    # Update workers table for canceled job with lifeLeft=-2
    # Note: the endDate MUST BE updated with now() so that it is not counted as alive workers
    logger.debug(JTM_SQL["update_lifeleft_enddate_workers_by_poolname"]
                 % dict(poolname=taskPool,
                        now=time.strftime("%Y-%m-%d %H:%M:%S")))
    db.execute(JTM_SQL["update_lifeleft_enddate_workers_by_poolname"]
               % dict(poolname=taskPool,
                      now=time.strftime("%Y-%m-%d %H:%M:%S")))
    db.close()
    logger.info("Pool, %s removed!" % taskPool)
    send_msg_callback(ch, method, props, 1, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)


# -------------------------------------------------------------------------------
def on_task_request(ch, method, props, body):
    """
    Event handler for processing request from jtm-* CLI tools
    """
    # Uncompress msg
    msgUnzipped = json.loads(zloads(body))

    logger.info("New task: {}".format(msgUnzipped))
    assert "task_type" in msgUnzipped, "Critical: need a task type!"
    taskType = msgUnzipped["task_type"]
    assert taskType in TASK_TYPE, "Critical: invalid task type: {}".format(msgUnzipped)
    innerTaskReqQ = JTM_INNER_REQUEST_Q

    # If pool is set, the tasks which use the queue name (=pool name)
    # will only be sent to the pool of workers
    # NOTE: poolName might be None even if "pool" in msgUnzipped
    if "pool" in msgUnzipped:
        poolJson = json.loads(json.dumps(msgUnzipped["pool"]))
        poolName = poolJson["name"] if "name" in poolJson and poolJson["name"] else None
        if poolName:
            # _jtm_inner_request_queue.<cluster_name>.jtm.<pool_name>
            innerTaskReqQ = JTM_INNER_REQUEST_Q + "." + poolName
        else:
            innerTaskReqQ = JTM_INNER_REQUEST_Q + ".small"

    # This allows jtm to keep task requests without any alive worker
    ch.exchange_declare(exchange=JTM_INNER_MAIN_EXCH,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)

    # Get queue length
    # taskQueueLen = ch.queue_declare(queue=innerTaskReqQ,
    #                                 durable=True,
    #                                 exclusive=False,
    #                                 auto_delete=True).method.message_count
    # logger.debug("#tasks queued = %d", taskQueueLen)

    ch.queue_declare(queue=innerTaskReqQ,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)
    ch.queue_bind(exchange=JTM_INNER_MAIN_EXCH,
                  queue=innerTaskReqQ,
                  routing_key=innerTaskReqQ)

    if taskType == "task":
        # TODO: innerTaskReqQ can be customized for different types of workers ==> done
        # ch.queue_bind() needs to be called with new queue name
        process_task_request(ch, method, props, msgUnzipped, innerTaskReqQ)

    elif taskType == "status":
        process_task_status(ch, method, props, int(msgUnzipped["task_id"]))

    elif taskType == "resource":
        # Note: this just returns the resource log file name. The file content should be converted
        #  into JSON format on the user side b/c the file might only be accessible from the user's
        #  machine.
        process_resource_log(ch, method, props, int(msgUnzipped["task_id"]))

    elif taskType == "kill":
        # TODO: need to check msgUnzipped["jtm_host_name"] to determine a way to kill a task per
        #  cluster/cloud
        process_task_kill(ch, method, props, msgUnzipped)

    elif taskType == "check_manager":
        # Just reply back to jtm-check-manager with '88'
        send_msg_callback(ch, method, props, 88, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)

    elif taskType == "check_worker":
        process_check_worker(ch, method, props, msgUnzipped)

    elif taskType == "remove_pool":
        # TODO: need to check msgUnzipped["jtm_host_name"] to determine a way to remove pool per
        #  cluster/cloud
        process_remove_pool(ch, method, props, msgUnzipped)

    else:
        logger.critical("Undefined task type: {}".format(taskType))
        send_msg_callback(ch, method, props, -1, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)
        logger.debug("sent reply: %s", str(-1))


# -------------------------------------------------------------------------------
def kill_all_workers(pool):
    """
    Send termination signal to all the workers
    :param pool: worker pool name (= taskqueue name) to kill
    :return:
    """
    # Remote broker (rmq.nersc.gov) connection open
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()
    ch.exchange_declare(exchange=JTM_INNER_MAIN_EXCH,
                        exchange_type="fanout",
                        durable=False,
                        auto_delete=False)

    msgContainer = {}
    msgContainer["task_type"] = TASK_TYPE["term"]
    msgContainer["task_queue"] = pool
    msgZipped = zdumps(json.dumps(msgContainer))

    assert pool.endswith(CNAME)
    try:
        ch.basic_publish(exchange=JTM_INNER_MAIN_EXCH,
                         routing_key=pool,
                         body=msgZipped)
    except Exception as detail:
        logger.exception("Exception: Failed to submit a TERM signal to the workers.")
        logger.exception("Detail: %s", str(detail))
        sys.exit(1)

    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def task_kill_thread():
    """
    Check "runs" table if any row with cancelled request (cancelled=1).
    Check if cancelled=1 and wid!=None and status!=-4(which means it's in running status)
    Then, send a kill request to the worker with task ID.

    """
    while 1:
        db = DbSqlMy(db=MYSQL_DB)
        # Get a list of task ids where cancelled -> requested but status != terminated
        tids = [int(i) for i in db.selectAs1Col(JTM_SQL["select_tids_runs_by_cancelled_and_wid"])]
        db.close()

        for tid in tids:
            try:
                # Get the wid and child pid
                db = DbSqlMy(db=MYSQL_DB)
                wid = db.selectAs1Col(JTM_SQL["select_workerid_workers_by_tid"] % dict(tid=tid))
                childPid = db.selectScalar(JTM_SQL["select_chilpid_runs_by_tid"] % dict(tid=tid))
                logger.debug("SQL: select_chilpid_runs_by_tid = %d" % int(childPid) if childPid else 0)
                logger.debug("SQL: select_workerid_workers_by_tid %d = %s" % (tid, wid))
                db.close()

                assert childPid > 0
                assert wid is not None

                # Send task id and process id to worker id
                send_task_kill_request(tid, wid[0], childPid)

                # Update runs table with "terminated" status
                db = DbSqlMy(db=MYSQL_DB)
                db.execute(JTM_SQL["update_runs_status_to_terminated_by_tid"]
                           % dict(tid=tid,
                                  sid=TASK_STATUS["terminated"],
                                  now=time.strftime("%Y-%m-%d %H:%M:%S")))

                # Note: do we have to update tasks with status = "terminated"?
                #  HB recv will check if it's terminated and update tasks table properly!
                #
                db.commit()
                db.close()
            except:
                logger.exception("Something goes wrong in task_kill_thread().")
                raise

        time.sleep(TASK_KILL_INTERVAL)


# -------------------------------------------------------------------------------
def zombie_worker_cleanup_thead():
    """
    Try to find cancelled slurm job and update workers table

    """
    while 1:
        db = DbSqlMy(db=MYSQL_DB)
        slurmJids = [int(i) for i in db.selectAs1Col(JTM_SQL["select_slurmjid_workers_by_lifeleft"])]
        db.close()
        if len(slurmJids) > 0:
            logger.info("Worker checking for %s" % str(slurmJids))
        for j in slurmJids:
            # Note: slurm dependent code!
            cmd = "sacct -j %d" % j
            so, se, ec = run_sh_command(cmd, live=True, log=logger)
            if ec == 0 and so.find("CANCELLED") != -1:
                db = DbSqlMy(db=MYSQL_DB)
                db.execute(JTM_SQL["update_workers_lifeleft_by_slurmjid"]
                           % dict(slurmjid=j,
                                  now=time.strftime("%Y-%m-%d %H:%M:%S")))
                db.commit()
                db.close()
            time.sleep(1)

        time.sleep(WORKER_KILL_INTERVAL)


