#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)
# pylint: disable=C0111,C0103,R0205

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
                      task request queue and task result queue in run time; Added jtm_host_name to
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
                      ex) task_id=228 --> 00/00/02;

    05.10.2019 5.2.0: Created single db connection for worker's hb recv;
    05.15.2019 5.2.1: Improved jtm-kill by sending kill msg to workers only from task_kill_proc();
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
                      ** single db conn open/close in recv_hb_from_worker_proc;
    06.13.2019 5.4.2: Multiprocessing -> threading;
                      recv_result -> multithreaded (default: 10)
    06.25.2019 5.4.4: Revered back to multiprocessing
                      (ref. https://stackoverflow.com/questions/3044580/multiprocessing-vs-threading-python)
    08.13.2019 5.4.5: Removed threads for on_result(); pika upgraded to v1.1.0;
    09.19.2019 5.4.6: Found missing db connection; Open and close db multiple times in recv_hb_from_worker_proc;
                      Delete tag 5.5.0 (= thread test)
                      Released as a stable production version;
    10.03.2019 5.4.7: Testing sending hb rates (worker sending: 1sec, manager recving: 3sec)

    11.21.2019 5.6.0: Set heartbeat=0 to prevent possible lost connection in BlockingConnection (in pika 0.9,
                      it was set to 580sec. In pika 1.0, it is set to 60sec)
    11.26.2019 5.6.1: Fixed bug for wrong number of workers -> needed to set nWorkersPerNode in workers table to 1;
                      Fixed process_check_worker() to get the correct number of workers alive
                      --> Added life_left>-0
                      to select_sum_nwpn_workers_by_jtm_host_name_enddate
                         select_sum_nwpn_workers_by_jtm_host_name_enddate
                         select_sum_nwpn_workers_by_poolname_enddate SQLs;
    12.09.2019 5.6.2: Updated to use mysql.connector;
    02.03.2020 5.6.3: Jtm-status fix; jtm-kill updates runs table for cancelled=1 but jtm-status still
                      checked only "status" field. So changed jtm-status to check "cancelled" field
                      to determine the status.
    02.12.2020 5.6.4: Updated log file permissions;
    02.14.2020 5.6.5: Updated to support custom charging account for knl;

    03.03.2020 5.6.7: Updated to use custom queue name instead of random in jtm-* interface;

"""
import multiprocessing as mp
import time
import json
from math import ceil
import uuid
import datetime
import sys
import pika
import shortuuid
import os
# import getpass
# import signal

from jaws_jtm.common import setup_custom_logger, logger
from jaws_jtm.lib.sqlstmt import JTM_SQL
from jaws_jtm.lib.rabbitmqconnection import RmqConnectionHB, send_msg_callback
from jaws_jtm.lib.dbutils import DbSqlMy
from jaws_jtm.lib.run import pad_string_path, make_dir, run_sh_command
from jaws_jtm.lib.msgcompress import zdumps, zloads

from jaws_jtm.config import JtmConfig
config = JtmConfig()
TASK_STATUS = config.constants.TASK_STATUS
VERSION = config.constants.VERSION
TASK_TYPE = config.constants.TASK_TYPE
WORKER_TYPE = config.constants.WORKER_TYPE
DONE_FLAGS = config.constants.DONE_FLAGS
HB_MSG = config.constants.HB_MSG
NUM_MANAGER_PROCS = config.constants.NUM_MANAGER_PROCS
PARENT_PROCESS_ID = os.getpid()  # parent process id

MYSQL_HOST = config.configparser.get("MYSQL", "host")
MYSQL_USER = config.configparser.get("MYSQL", "user")
MYSQL_PORT = config.configparser.getint("MYSQL", "port")
MYSQL_DB = config.configparser.get("MYSQL", "db")
JTM_INNER_MAIN_EXCH = config.configparser.get("JTM", "jtm_inner_main_exch")
CNAME = config.configparser.get("SITE", "instance_name")
JTM_TASK_RESULT_Q = config.configparser.get("JTM", "jtm_task_result_q")
JTM_TASK_REQUEST_Q = config.configparser.get("JTM", "jtm_task_request_q")
JGI_JTM_MAIN_EXCH = config.configparser.get("JTM", "jgi_jtm_main_exch")
JTM_INNER_REQUEST_Q = config.configparser.get("JTM", "jtm_inner_request_q")
JTM_INNER_RESULT_Q = config.configparser.get("JTM", "jtm_inner_result_q")
WORKER_HB_RECV_INTERVAL = config.configparser.getfloat("JTM", "worker_hb_recv_interval")
JTM_WORKER_POISON_EXCH = config.configparser.get("JTM", "jtm_worker_poison_exch")
JTM_WORKER_POISON_Q = config.configparser.get("JTM", "jtm_worker_poison_q")
JTM_TASK_KILL_EXCH = config.configparser.get("JTM", "jtm_task_kill_exch")
JTM_TASK_KILL_Q = config.configparser.get("JTM", "jtm_task_kill_q")
TASK_STAT_UPDATE_INTERVAL = config.configparser.getfloat("JTM", "task_stat_update_interval")
JTM_CLIENT_HB_EXCH = config.configparser.get("JTM", "jtm_client_hb_exch")
JTM_WORKER_HB_EXCH = config.configparser.get("JTM", "jtm_worker_hb_exch")
WORKER_KILL_INTERVAL = config.configparser.getfloat("JTM", "worker_kill_interval")
TASK_KILL_INTERVAL = config.configparser.getfloat("JTM", "task_kill_interval")
CLIENT_HB_SEND_INTERVAL = config.configparser.getfloat("JTM", "client_hb_send_interval")
NUM_RESULT_RECV_THREADS = config.configparser.getint("JTM", "num_result_recv_threads")
NUM_PROCS_CHECK_INTERVAL = config.configparser.getfloat("JTM", "num_procs_check_interval")
ENV_ACTIVATION = config.configparser.get("JTM", "env_activation")
RESULT_RECEIVE_INTERVAL = config.configparser.getfloat("JTM", "result_receive_interval")
RUNS_INFO_UPDATE_WAIT = config.configparser.getfloat("JTM", "runs_info_update_wait")
WORKER_INFO_UPDATE_WAIT = config.configparser.getfloat("JTM", "worker_info_update_wait")


# --------------------------------------------------------------------------------------------------
# Globals
# --------------------------------------------------------------------------------------------------
NUM_TOTAL_WORKERS = mp.Value('i', 0)


# --------------------------------------------------------------------------------------------------
def recv_hb_from_worker_proc(hb_queue_name, log_dest_dir, b_resource_log):
    """
    Receive hearbeats from all the workers running
    :param hb_queue_name: hb queue name
    :param hb_send_proc_handle: hb sending process handle
    :param log_dest_dir: log destination path
    :param b_resource_log: if true, create resource log file
    :return:
    """

    # Todo: change to event driven by ch.basic_consume and ch.start_consuming
    # Remote broker (mq.nersc.gov)
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()

    exch_name = JTM_WORKER_HB_EXCH

    ch.exchange_declare(exchange=exch_name,
                        exchange_type="direct",
                        durable=False,
                        auto_delete=False)

    # This queue can be declared from workers first.
    try:
        ch.queue_declare(queue=hb_queue_name,
                         durable=False,
                         exclusive=False,
                         auto_delete=True)
    except Exception as detail:
        logger.exception("Exception: The queue, %s is in use.", hb_queue_name)
        logger.exception("Detail: %s", str(detail))
        # ch.stop_consuming()
        # ch.close()
        # conn.close()
        # sys.exit(1)
        raise

    # Queue binding for the queue to receive hb from workers
    ch.queue_bind(exchange=exch_name,
                  queue=hb_queue_name,
                  routing_key=hb_queue_name)

    b_is_msg_cleared = False  # all stacked messages are processed or not
    b_is_worker_found = False  # is any alive worker
    max_worker_check_count = 0  # max number of checking workers
    interval = config.configparser.getfloat("JTM", "client_hb_recv_interval")

    # Todo: need to open/close in place to evade from locking. Any other solutions?
    # Note: tried connection pool but not working as expected
    # db = DbSqlLite(dbpath="jtm.db")
    # db = DbSqlMy(db=MYSQL_DB)

    while 1:
        worker_ids_dict = {}
        try:
            method_frame, header_frame, body = ch.basic_get(queue=hb_queue_name, auto_ack=True)

            # Get all the hb messages from workers
            # if method_frame:
            #    for i in range(int(method_frame.message_count)):
            #        method_frame, header_frame, body = ch.basic_get(queue=hb_queue_name, auto_ack=True)
        except Exception as detail:
            logger.exception("Exception: Failed to get a message from %s.", hb_queue_name)
            logger.exception("Detail: %s", str(detail))
            # ch.stop_consuming()
            # ch.close()
            # conn.close()
            # os._exit(1)
            raise

        if body and not b_is_msg_cleared:
            # To clear up any heartbeat messages left in the heartbeat queue.
            for i in range(int(method_frame.message_count)):
                _, _, _ = ch.basic_get(queue=hb_queue_name, auto_ack=True)

            b_is_msg_cleared = True

        elif body and b_is_msg_cleared:
            # To deal with the heartbeats from the worker(s). The workers send
            # each hostname and pid. Using "dict", get the number of unique
            # pids of the running workers.
            msg_unzipped = json.loads(zloads(body))
            # type conversion and sort by key
            msg_unzipped = {int(k): v for k, v in msg_unzipped.items()}

            # worker id is used to collect unique root_proc_id (= num of workers)
            a_worker_id = msg_unzipped[HB_MSG["worker_id"]]
            worker_ids_dict[a_worker_id] = msg_unzipped

            # NOTE: Workers send it"s hb interval to the client in the msg packet.
            # Set the checking heartbeat as (jtm-worker"s sending heartbeat
            # interval value * intervalIncRate)

            # Worker's hb should be more than one so consume all the hb's
            for i in range(int(method_frame.message_count)):
                method_frame, header_frame, body = ch.basic_get(queue=hb_queue_name, auto_ack=True)

                msg_unzipped = json.loads(zloads(body))
                msg_unzipped = {int(k): v for k, v in msg_unzipped.items()}
                a_worker_id = msg_unzipped[HB_MSG["worker_id"]]
                worker_ids_dict[a_worker_id] = msg_unzipped

            # for k, v in worker_ids_dict.iteritems():
            for k, v in worker_ids_dict.items():  # py3
                task_id = v[HB_MSG["task_id"]]
                root_proc_id = v[HB_MSG["root_pid"]]
                child_proc_id = v[HB_MSG["child_pid"]]
                a_worker_id = v[HB_MSG["worker_id"]]
                slurm_job_id = v[HB_MSG["slurm_jobid"]]
                worker_type = v[HB_MSG["worker_type"]]
                end_datetime = v[HB_MSG["end_date"]]
                life_left = v[HB_MSG["life_left"]]
                mem_per_node = v[HB_MSG["mem_per_node"]]
                mem_per_core = v[HB_MSG["mem_per_core"]]
                num_cores = v[HB_MSG["num_cores"]]
                job_time = v[HB_MSG["job_time"]]
                clone_time = v[HB_MSG["clone_time_rate"]]
                host_name = v[HB_MSG["host_name"]]
                jtm_host_name = v[HB_MSG["jtm_host_name"]]
                ip_addr = v[HB_MSG["ip_address"]]
                pool_name = v[HB_MSG["pool_name"]]
                # num_worker_on_the_node = v[HB_MSG["nwpn"]]  # Todo: now with nwpn, num_worker_on_the_node = 1 always

                if b_resource_log:
                    logger.resource(v)

                if pool_name:
                    queue_name = JTM_INNER_REQUEST_Q + "." + pool_name
                else:
                    queue_name = ""

                # Fixme: bytearray index out of range EXCEPTION with gpdb23
                success = False

                # Todo: still need to do "while" for checking db connection?
                while success is not True:
                    success = True
                    try:
                        db = DbSqlMy(db=MYSQL_DB)
                        # This increases AUTO_INCREMENT field even there is no insertion
                        # db.execute("""INSERT IGNORE INTO workers (a_worker_id, slurm_job_id, worker_type)
                        #     VALUES ("%(worker_id)s", %(slurm_jid)d, %(worker_type)d)
                        # """ % dict(worker_id=a_worker_id,
                        #            slurm_jid=slurm_job_id,
                        #            worker_type=worker_type))
                        bExists = db.selectScalar(JTM_SQL["select_exists_workers_by_workerid"]
                                                  % dict(worker_id=a_worker_id))

                        if not bExists:
                            # Todo: remove num_worker_on_the_node insertion
                            db.execute(JTM_SQL["insert_workers_workerid_slurmjobid"]
                                       % dict(worker_id=a_worker_id,
                                              slurm_jid=slurm_job_id,
                                              worker_type=worker_type,
                                              # nwpn=num_worker_on_the_node,
                                              nwpn=1))
                        else:
                            # Todo: stress test needed! Test if nworkers > 1000
                            # if worker_type == WORKER_TYPE["manual"]
                            #     db.execute(JTM_SQL["update_workers_enddate_lifeleft_by_workerid"]
                            #                % dict(worker_id=a_worker_id, now=end_datetime))
                            # elif worker_type == WORKER_TYPE["static"]:
                            #     pass
                            # elif worker_type == WORKER_TYPE["dynamic"]:
                            #     pass
                            # Todo: can only update end_datetime and life_left after 1st insert
                            db.execute(JTM_SQL["update_workers_enddate_lifeleft_by_workerid"]
                                       % dict(worker_id=a_worker_id,
                                              now=end_datetime,
                                              life_left=life_left,
                                              mpn=mem_per_node if worker_type != 0 else 0,
                                              mpc=mem_per_core if worker_type != 0 else 0,
                                              num_cores=num_cores if worker_type != 0 else 0,
                                              job_time=job_time,
                                              clone_rate=clone_time if worker_type != 0 else 0,
                                              host_name=host_name,
                                              jtm_host_name=jtm_host_name,
                                              ipaddr=ip_addr,
                                              pool_name=queue_name,
                                              slurm_jid=slurm_job_id))

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

                if task_id > 0 and root_proc_id != child_proc_id:  # if there is a user process started
                    datetime_str = datetime.datetime.now().strftime("%Y-%m-%d")
                    if log_dest_dir:
                        log_dir_name = os.path.join(log_dest_dir, "resource")
                    else:
                        log_dir_name = "%s/resource" % (os.getcwd())

                    padded_dir_str = pad_string_path(task_id, depth=3)  # 226 --> 00/00/02
                    make_dir(os.path.join(log_dir_name, padded_dir_str))
                    resourceLogFileName = "%s/%s/jtm_resource_%d_%s.log" \
                                          % (log_dir_name, padded_dir_str, task_id, datetime_str)
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
                    os.chmod(resourceLogFileName, 0o777)

                    # Update runs table with resource usage for a task
                    # rsc = ",".join([str(i) for i in collections.OrderedDict(sorted(v.items())).values()])
                    # Update runs table
                    # db.execute(JTM_SQL["update_runs_resources_by_tid"]
                    #            % dict(rsc=rsc, task_id=task_id))

                    try:
                        # Update tasks table with "running" status == 2 if status is still 0 or 1
                        db = DbSqlMy(db=MYSQL_DB)
                        # db.execute(JTM_SQL["update_runs_status_workerid_by_tid"]
                        #            % dict(status_id=TASK_STATUS["running"],
                        #                   task_id=task_id,
                        #                   worker_id=a_worker_id,
                        #                   child_pid=child_proc_id))
                        # db.execute(JTM_SQL["insert_ignore_workers_by_wid_slurmjid"]
                        #            % dict(worker_id=a_worker_id,
                        #                   slurm_jid=slurm_job_id))

                        # Update runs table for a task
                        # Todo: really need to store full path to the log?
                        # logger.debug(db.selectAll("select * from runs where task_id=%d" % task_id))
                        db.execute(JTM_SQL["update_runs_status_by_taskid"]
                                   % dict(status_id=TASK_STATUS["running"],
                                          task_id=task_id,
                                          worker_id=a_worker_id,
                                          child_pid=child_proc_id,
                                          resource_log=resourceLogFileName))
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

                job_time_to_minute = 0
                if job_time:
                    job_time_to_minute = hms_to_m(job_time)

                rate = float(life_left) / job_time_to_minute if life_left > 0 else 0.0

                if worker_type == 1:  # only for static workers
                    logger.debug("worker type {}, wid {}, jobtime {}, jobtime in minute {}, "
                                 "Life left in minute {}, computedRate {}, cloneTimeRate {}".format(worker_type,
                                                                                                    a_worker_id,
                                                                                                    job_time,
                                                                                                    job_time_to_minute,
                                                                                                    life_left,
                                                                                                    rate,
                                                                                                    clone_time))

                # Static workers cloning is determined by the lifeleft
                # if job_time and worker_type > 0 and slurm_job_id > 0 and rate <= float(clone_time):
                if job_time and worker_type == 1 and slurm_job_id > 1 and rate <= float(clone_time):
                    db = DbSqlMy(db=MYSQL_DB)
                    cloneCnt = db.selectScalar(JTM_SQL["select_clonecnt_workers_by_workerid"]
                                               % dict(worker_id=a_worker_id))
                    db.close()

                    logger.debug("Clone count = %s" % cloneCnt)

                    if int(cloneCnt) == 0:
                        try:
                            db = DbSqlMy(db=MYSQL_DB)
                            # db.execute(JTM_SQL["update_workers_clonecnt_by_workerid"]
                            #            % dict(worker_id=a_worker_id))

                            # To handle multiple static workers on a same node (-nw option)
                            # increase clonecnt by job id 01092019
                            db.execute(JTM_SQL["update_workers_clonecnt_by_slurmjobid"]
                                       % dict(slurmjobid=slurm_job_id))
                            db.commit()
                            db.close()
                        except Exception as e:
                            logger.critical(e)
                            logger.critical("Failed to update workers table for clonecnt.")
                            raise

                        logger.info("Send cloning signal to the worker, {}".format(a_worker_id))

                        # Todo: if multiple workers are in a node, need to send this reproduce command
                        #  to only one of those workers
                        #  only static worker spawn itself automatically
                        if worker_type == WORKER_TYPE["static"]:
                            send_reproduce_or_die(a_worker_id, 0, 1)

                # Dynamic workers cloning and terminating are determined by
                # total # tasks in the queue - number_workers_alive + workers_requested
                #
                # 11.14.2018
                # no automatic cloning of dynamic workers for now
                # user specified custom or base pool will be created and the pool will be disappeared after jobtime
                #
                # if worker_type == WORKER_TYPE["dynamic"]:
                #     # Todo: need num clones option, default: 2
                #     # Todo: need to take #tasks in the queue into consideration
                #     send_reproduce_or_die(a_worker_id, 0, 2)
                ########################################################################################################

            #
            # Set unresponsive workers as dead
            #
            # Collect worker_ids from hb
            alive_worker_id_list = []
            # for k, v in worker_ids_dict.iteritems():
            for k, v in worker_ids_dict.items():  # py3
                alive_worker_id_list.append(v[HB_MSG["worker_id"]])
            # logger.debug("# unique worker IDs in HB = %d" % len(alive_worker_id_list))

            # Collect worker_id which are still set as alive from workers table
            success = False
            while success is not True:
                success = True
                try:
                    db = DbSqlMy(db=MYSQL_DB)
                    # live_worker_id_list = db.selectAll(JTM_SQL["select_workerid_workers_by_lifeleft"]
                    live_worker_id_list = db.selectAll(JTM_SQL["select_workerid_workers_by_lifeleft_jtmhostname"]
                                                       % dict(jtm_host_name=jtm_host_name))
                    live_worker_id_list = [i[0] for i in live_worker_id_list]
                    # logger.debug("# of unique worker IDs in table = %d" % len(live_worker_id_list))

                    # Check & update workers table
                    # set life_left as -1 for dead workers
                    for w in live_worker_id_list:
                        if w not in alive_worker_id_list:
                            db.execute(JTM_SQL["update_workers_lifeleft_by_workerid"]
                                       % dict(worker_id=w, jtm_host_name=jtm_host_name))

                            # Todo: should delete hb queue?

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

            # NUM_TOTAL_WORKERS.value = len(alive_worker_id_list)
            db = DbSqlMy(db=MYSQL_DB)
            # alive_total_num_workers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_lifeleft"])
            alive_total_num_workers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_lifeleftt_jtmhostname"]
                                                      % dict(jtm_host_name=jtm_host_name))
            # logger.debug(JTM_SQL["select_sum_nwpn_workers_by_lifeleftt_jtmhostname"]
            #              % dict(jtm_host_name=jtm_host_name))

            db.close()
            NUM_TOTAL_WORKERS.value = int(alive_total_num_workers) if alive_total_num_workers else 0
            logger.debug("# workers: in hb=%d, in table=%d, alive+requested=NUM_TOTAL_WORKERS=%d"
                         % (len(alive_worker_id_list), len(live_worker_id_list), NUM_TOTAL_WORKERS.value))

            if NUM_TOTAL_WORKERS.value > 0:
                b_is_worker_found = True
                max_worker_check_count = 0  # reinitialize

        elif not body and b_is_msg_cleared and b_is_worker_found:
            # If we are here, unfortunately, we lost all the workers that we've been using.
            max_worker_check_count += 1
            NUM_TOTAL_WORKERS.value = 0
            logger.info("Waiting for worker(s)...")

            # NOTE: 2013.09.05 To prevent from failing to detect workers
            # Increase interval. default=1.2
            # Todo: still need this? ==> 11.13.2018 removed
            # intervalIncRate = intervalIncRate * CLIENT_HB_RECEIVE_INT_INC_RATE
            WORKER_HB_CHECK_MAX_COUNT = config.configparser.getint("JTM", "worker_hb_check_max_count")
            if WORKER_HB_CHECK_MAX_COUNT != 0 and \
                    max_worker_check_count > WORKER_HB_CHECK_MAX_COUNT:  # hit the max checking limit
                # Close connection and kill parent and itself
                raise OSError(2, 'Worker not found')

            # If there no workers alive after 5 checks, set life_left to -1 for all
            if max_worker_check_count == 3:
                try:
                    db = DbSqlMy(db=MYSQL_DB)
                    db.execute(JTM_SQL["update_workers_lifeleft_for_last"])
                    db.commit()
                    db.close()
                except Exception as e:
                    logger.critical(e)
                    logger.critical("Failed to update workers table for lifeleft.")
                    raise

        # Todo: Need to start with shorter interval like (0.5-1 sec) for about 30 sec
        #  from the beginning and then use the user specified interval
        #  ==> 11.13.2018 just wait CLIENT_HB_RECV_INTERVAL
        # time.sleep(interval)

        # To fix connection lost
        # Ref) https://github.com/pika/pika/issues/1224
        try:
            conn.process_data_events(time_limit=float(interval))
        except Exception as e:
            logger.critical(e)
            logger.critical("RMQ connection lost.")
            # ch.stop_consuming()
            # ch.close()
            # conn.close()
            # os._exit(1)
            raise

    # unreachable
    # db.close()
    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def send_hb_to_worker_proc():
    """
    Broadcast heartbeat to all workers
    Ref) http://www.rabbitmq.com/tutorials/tutorial-three-python.html
    """
    # Remote broker (mq.nersc.gov)
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()
    exch_name = JTM_CLIENT_HB_EXCH

    ch.exchange_declare(exchange=exch_name,
                        exchange_type="topic",
                        durable=False,
                        auto_delete=False)

    msg_to_send_dict = {}
    msg_to_send_dict["task_type"] = TASK_TYPE["hb"]
    msg_zipped = zdumps(json.dumps(msg_to_send_dict))

    try:
        while 1:
            ch.basic_publish(exchange=exch_name,
                             routing_key="*." + CNAME,  # all workers with CNAME can hear it
                             body=msg_zipped)
            conn.process_data_events(time_limit=CLIENT_HB_SEND_INTERVAL)
    except Exception as e:
        logger.critical("Something wrong in send_hb_to_worker_proc(): %s", e)
        raise

    # unreachable
    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def recv_result_from_workers_proc():
    """
    Receive hearbeats from all the workers running
    """
    # Remote broker (mq.nersc.gov)
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()

    exch_name = JTM_INNER_MAIN_EXCH
    inner_result_queue_name = JTM_INNER_RESULT_Q

    # Default; exch_name = jgi_jtm_inner_main_exchange
    ch.exchange_declare(exchange=exch_name,
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
        ch.queue_declare(queue=inner_result_queue_name,
                         durable=True,
                         exclusive=False,
                         auto_delete=True)

    except Exception as detail:
        logger.exception("Exception: The queue, %s is in use.", JTM_TASK_RESULT_Q)
        logger.exception("Detail: %s", str(detail))
        # ch.stop_consuming()
        # ch.close()
        # conn.close()
        # # sys.exit(1)
        # os._exit(1)
        raise

    # Queue binding for the queue to receive hb from workers
    ch.queue_bind(exchange=exch_name,
                  queue=inner_result_queue_name,
                  routing_key=inner_result_queue_name)

    # Todo: change to a threaded version
    ch.basic_qos(prefetch_count=NUM_RESULT_RECV_THREADS)

    # NOTE: the below methods might cause error in consuming messages which are already queued
    # ch.basic_consume(lambda ch, method, properties, body: recv_result_on_result(ch, method, properties, body, dbpool),
    #                  queue=inner_result_queue_name,
    #                  auto_ack=False)
    # recv_result_on_result_callback = functools.partial(recv_result_on_result, args=(dbpool))
    # ch.basic_consume(recv_result_on_result_callback, inner_result_queue_name, auto_ack=False)

    # OLD
    try:
        ch.basic_consume(queue=inner_result_queue_name,
                         on_message_callback=recv_result_on_result)
    except Exception as detail:
        logger.exception("Exception: The queue, %s is in use.", JTM_TASK_RESULT_Q)
        logger.exception("Detail: %s", str(detail))
        # ch.stop_consuming()
        # ch.close()
        # conn.close()
        # # sys.exit(1)
        # os._exit(1)
        raise

    # NEW
    # threads = []
    # on_result_callback = functools.partial(on_result, args=(conn, threads))
    # ch.basic_consume(inner_result_queue_name, on_result_callback)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        ch.stop_consuming()
        raise

    # Wait for all to complete
    # for thread in threads:
    #     thread.join()

    # Fixme: even there are a result in _jtm_inner_result_queue.* queue, it is not consumed.
    # while 1:
    #     method_frame, header_frame, body = ch.basic_get(queue=inner_result_queue_name, auto_ack=False)
    #     if body:
    #         msg_unzipped = json.loads(zloads(body))
    #         done_flag = int(msg_unzipped["done_flag"])
    #         task_id = int(msg_unzipped["task_id"])
    #         ret_msg = msg_unzipped["ret_msg"]
    #         a_worker_id = msg_unzipped["worker_id"]
    #         host_name = msg_unzipped["host_name"]
    #
    #         logger.debug("Result received: {}".format(msg_unzipped))

    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def nack_message(ch, delivery_tag):
    """Note that `ch` must be the same pika channel instance via which
    the message being ACKed was retrieved (AMQP protocol constraint).
    """
    if ch.is_open:
        # ch.basic_ack(delivery_tag)
        ch.basic_reject(delivery_tag=delivery_tag, requeue=False)
    else:
        # Channel is already closed, so we can't ACK this message;
        # log and/or do something that makes sense for your app in this case.
        pass


# -------------------------------------------------------------------------------
def ack_message(ch, delivery_tag):
    """Note that `ch` must be the same pika channel instance via which
    the message being ACKed was retrieved (AMQP protocol constraint).
    """
    if ch.is_open:
        ch.basic_ack(delivery_tag)
    else:
        # Channel is already closed, so we can't ACK this message;
        # log and/or do something that makes sense for your app in this case.
        pass


# -------------------------------------------------------------------------------
def recv_result_on_result(ch, method, props, body):
    """
    recv_result_from_workers_proc's basic_consume callback
    :param ch: channel
    :param method: method_frame
    :param props: pika connection property
    :param body: message received
    :return:
    """
    if body:
        # try:
        #     msg_unzipped = json.loads(zloads(body))
        #     done_flag = int(msg_unzipped["done_flag"])
        #     task_id = int(msg_unzipped["task_id"])
        #     ret_msg = msg_unzipped["ret_msg"]
        #     a_worker_id = msg_unzipped["worker_id"]
        #     host_name = msg_unzipped["host_name"]
        # except KeyError as ke:
        #     logger.exception("Key error: {}".format(ke))
        #     logger.info("msg = {}".format(msg_unzipped))
        #     ch.stop_consuming()
        #     ch.close()
        #     os._exit(1)
        # except Exception as e:
        #     raise
        try:
            msg_unzipped = json.loads(zloads(body))
        except Exception as e:
            logger.critical(e)
            logger.exception("Failed to load returned result message, {}".format(body))
            # ch.stop_consuming()
            # ch.close()
            # os._exit(1)
            raise

        assert "task_id" in msg_unzipped
        task_id = int(msg_unzipped["task_id"])
        done_flag = int(msg_unzipped["done_flag"]) if "done_flag" in msg_unzipped else -2
        ret_msg = msg_unzipped["ret_msg"] if "ret_msg" in msg_unzipped else ""
        a_worker_id = msg_unzipped["worker_id"] if "worker_id" in msg_unzipped else ""
        host_name = msg_unzipped["host_name"] if "host_name" in msg_unzipped else ""

        # logger.debug("Result received: {} {}".format(method, props))
        logger.info("Result received: {}".format(msg_unzipped))

        if ret_msg != "hb":
            logger.debug("Update tasks {} with {}".format(task_id, done_flag))
            try:
                db = DbSqlMy(db=MYSQL_DB)
                db.execute(JTM_SQL["update_tasks_doneflag_by_taskid"]
                           % dict(task_id=task_id,
                                  done_flag=done_flag))
                db.commit()
                db.close()
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed to update tasks table for doneflag.")
                raise

            # Check the task status code
            if done_flag > 0:
                task_status_int = TASK_STATUS["success"]  # 4
            else:
                if done_flag == -1:
                    task_status_int = TASK_STATUS["outputerror"]
                elif done_flag == -2:
                    task_status_int = TASK_STATUS["failed"]
                elif done_flag == -3:
                    task_status_int = TASK_STATUS["outofresource"]
                elif done_flag == -4:
                    task_status_int = TASK_STATUS["terminated"]
                elif done_flag == -6:
                    task_status_int = TASK_STATUS["timeout"]
                else:
                    logger.critical("Unknown return code {}".format(done_flag))
                    # ch.stop_consuming()
                    # ch.close()
                    # sys.exit(1)
                    raise OSError(2)

            # This seems like resolving runs table lock issue
            # issue: status is not changed to 4 after the update
            # Todo: need to improve
            #############################################
            ch._connection.sleep(RESULT_RECEIVE_INTERVAL)
            #############################################

            logger.debug("Update runs for task id = {} with {} for worker {}".format(task_id,
                                                                                     task_status_int,
                                                                                     a_worker_id))

            # Sometimes workerId2 ==> 0
            # so wait a little bit if that happened
            a_worker_id_to_check = 0
            db = DbSqlMy(db=MYSQL_DB)
            rows = db.selectAll(JTM_SQL["select_workerid2_workers_by_wid"]
                                % dict(worker_id=a_worker_id))
            db.close()

            try:
                a_worker_id_to_check = int(rows[0][0])
            except Exception:
                a_worker_id_to_check = 0

            # Just in case
            while a_worker_id_to_check == 0:
                db = DbSqlMy(db=MYSQL_DB)
                rows = db.selectAll(JTM_SQL["select_workerid2_workers_by_wid"]
                                    % dict(worker_id=a_worker_id))
                db.close()
                try:
                    a_worker_id_to_check = int(rows[0][0])
                except Exception:
                    a_worker_id_to_check = 0

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
                               % dict(status_id=task_status_int,
                                      wid2=a_worker_id_to_check,
                                      now=time.strftime("%Y-%m-%d %H:%M:%S"),
                                      task_id=task_id))
                    logger.debug("status {}, worker id {}, taskid {}".format(task_status_int,
                                                                             a_worker_id_to_check,
                                                                             task_id))
                    db.commit()
                    # logger.debug(db.selectAll("select * from runs where taskId=%d" % task_id))
                    db.close()

                except Exception as e:
                    logger.critical(e)
                    logger.critical("Failed to update runs table for status and workerid2.")
                    logger.debug("Retry to update runs table for status and workerid2.")
                    # raise
                    success = False
                    logger.debug("update_runs_status_workerid2_by_taskid_2 sleep for %f" % RUNS_INFO_UPDATE_WAIT)
                    # time.sleep(RUNS_INFO_UPDATE_WAIT)
                    ch._connection.sleep(RUNS_INFO_UPDATE_WAIT)

            # Print report
            if done_flag == DONE_FLAGS["success"]:  # 1
                logger.info("Task %s --> Success on worker/host, %s/%s",
                            task_id, a_worker_id, host_name)

            elif done_flag == DONE_FLAGS["success with correct output file(s)"]:  # 2
                logger.info("Task %s --> Success with valid output(s) on worker/host, %s/%s",
                            task_id, a_worker_id, host_name)

            elif done_flag == DONE_FLAGS["failed to check output file(s)"]:  # -1
                logger.info("Task %s --> %s, worker/host: %s/%s",
                            task_id, ret_msg, a_worker_id, host_name)

            elif done_flag == DONE_FLAGS["failed to run user command"]:  # -2
                logger.info("Task %s --> Failed with non-zero exit code. stdout = %s, worker/host: %s/%s",
                            task_id, ret_msg, a_worker_id, host_name)

            elif done_flag == DONE_FLAGS["failed with out-of-mem"]:  # -3
                pass

            elif done_flag == DONE_FLAGS["failed with user termination"]:  # -4
                logger.info("Task %s --> Failed by user termination. stdout = %s, worker/host: %s/%s",
                            task_id, ret_msg, a_worker_id, host_name)

            elif done_flag == DONE_FLAGS["failed with timeout"]:  # -6
                logger.info("Task %s --> Failed with timeout. stdout = %s, worker/host: %s/%s",
                            task_id, ret_msg, a_worker_id, host_name)

            else:
                logger.warning("Cannot recognize the return code: %d" % done_flag)
                logger.info("Task %s --> worker/host: %s/%s", task_id, a_worker_id, host_name)
                # sys.exit(1)
                raise OSError(2, 'Cannot recognize the DONE_FLAG return code')

        else:
            # If it is not a result msg, return it back to the exchange
            ch.basic_reject(delivery_tag=method.delivery_tag, requeue=False)

        # After this the result message will be deleted from RabbitMQ
        # If this worker crashes while running a user command, this task will
        # be sent to other workers available
        ch.basic_ack(delivery_tag=method.delivery_tag)

    else:
        logger.critical("No data found from the result message.")
        # sys.exit(1)
        raise OSError(2, 'No data found from the result message')


# -------------------------------------------------------------------------------
def process_task_request(ch, method, props, msg, inner_task_request_queue):
    """
    Get task request from jtm-submit and send it to a worker
    :param ch:
    :param method:
    :param props:
    :param msg: unzipped dict
    :param inner_task_request_queue: inner task queue name (jgi-task-manager --> worker)
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

    user_task_cmd = msg["command"]
    task_type = msg["task_type"]
    output_file = msg["output_files"] if "output_files" in msg else ""  # comma separated list ex) "a.out,b.out,c.out"
    output_dir = msg["output_dir"] if "output_dir" in msg else ""
    stdout_file = msg["stdout"] if "stdout" in msg else ""
    stderr_file = msg["stderr"] if "stderr" in msg else ""
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
    b_failed_to_request_worker = False

    if "pool" in msg and "name" in msg["pool"] and "time" in msg["pool"]:
        # {u'resource': u'cori', u'name': u'test', u'size': 1}
        # ==> {"resource": "cori", "name": "test", "size": 1}
        pool_spec_json_str = json.loads(json.dumps(msg["pool"]))

        # Set default values defined in Config.py
        assert "name" in pool_spec_json_str
        pool_name = pool_spec_json_str["name"]
        assert pool_name is not None
        pool_cluster = config.configparser.get("JTM", "cluster")
        pool_ncpus = config.configparser.getint("JTM", "ncpus")
        pool_mem = config.configparser.get("JTM", "mempernode")
        pool_constraint = config.configparser.get("SLURM", "constraint")
        pool_charge_account = config.configparser.get("SLURM", "charge_accnt")
        pool_qos = config.configparser.get("SLURM", "qos")

        # Note: pool size = num_nodes_to_request * num_workers_per_node
        num_workers_per_node = config.configparser.getint("JTM", "num_workers_per_node")
        num_nodes_to_request = config.configparser.getint("JTM", "nnodes")

        # Worker type is restricted to "dynamic" for now.
        # Todo: add the feature to remove the custom pool by user task to support "static"
        #   workers in custom pool
        # worker_type = "dynamic"

        # NOTE: if pool_size is set, which means user sets the number of workers needed,
        # n number of jtm dynamic workers will be created, and the pool of n dynamic
        # workers won't be increased over n. The workers will be terminated if there is
        # no tasks requested
        #
        # if not,
        # only one dynamic worker will be created. The worker will be terminated if there is
        # no tasks in the queue for a specified time duration

        # if "size" in pool_spec_json_str:  # pool size = the number workers in the pool
        #     pool_size = int(pool_spec_json_str["size"])
        if "cluster" in pool_spec_json_str:  # cluster/clouds name
            pool_cluster = pool_spec_json_str["cluster"]
        if "time" in pool_spec_json_str:  # wallclocktime request
            pool_time = pool_spec_json_str["time"]
        if "cpu" in pool_spec_json_str:  # number of cores request
            pool_ncpus = int(pool_spec_json_str["cpu"])
        if "mem" in pool_spec_json_str:  # node memory request
            pool_mem = pool_spec_json_str["mem"]
        if "constraint" in pool_spec_json_str:  # [haswell | knl | skylake]
            pool_constraint = pool_spec_json_str["constraint"]
        if "qos" in pool_spec_json_str:
            pool_qos = pool_spec_json_str["qos"]  # ["genepool_special", "genepool_shared", "jgi_shared", "jgi_exvivo"]
        if "account" in pool_spec_json_str:
            pool_charge_account = pool_spec_json_str["account"]  # for example, gtrqc for skylake, fungalp for the rest

        # Todo: This is not used for now.
        #     #  If shared=0 is set from jtm-submit and if the pool needs to be terminated forcefully,
        #     #  jtm-manager should send a poison to the worker(s
        # if "shared" in pool_spec_json_str:
        #     pool_shared = pool_spec_json_str["shared"]
        if "nwpn" in pool_spec_json_str:
            num_workers_per_node = int(pool_spec_json_str["nwpn"])
        if "node" in pool_spec_json_str:
            num_nodes_to_request = int(pool_spec_json_str["node"])

        assert (len(user_task_cmd) <= 1024)
        assert (len(output_file) <= 1024)

        # Create pool
        # Step
        # 1. check the number if workers in the custom pool
        # 2. if 0, create a pool
        #    else send tasks to the pool

        # Todo: maintain the number of workers sbatched -->
        #   num_worker_to_add = pool_size - num_live_worker_in_pool - nWorkerSbatched
        db = DbSqlMy(db=MYSQL_DB)
        # num_live_worker_in_pool = db.selectScalar(JTM_SQL["select_count_workers_by_poolname_enddate"]
        #                                     % dict(pool_name=inner_task_request_queue,
        #                                            hbinterval=WORKER_HB_RECV_INTERVAL * 2))

        # Note: This is also to check the number of sbatched workers in the case of scatter operation in Cromwell
        # If the # of sbatched dynamic worker(s) is less than the requested pool size,
        # then do sbatch for the rest (=num_worker_to_add)
        #
        # num_live_worker_in_pool = db.selectScalar(JTM_SQL["select_count_workers_by_poolname"]
        #                                     % dict(pool_name=inner_task_request_queue,
        #                                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))

        # TODO: slurm jid --> sacct --> doublecheck the status of node allocation and worker
        #  if no real alive workers --> request new node
        slurm_jid_list = db.selectAll(JTM_SQL["select_distinct_jid_workers_by_poolname"]
                                      % dict(pool_name=inner_task_request_queue,
                                             hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("slurm id for {}: {}".format(inner_task_request_queue, slurm_jid_list))

        num_slurm_jid_in_pool = db.selectScalar(JTM_SQL["select_count_distinct_jid_workers_by_poolname"]
                                                % dict(pool_name=inner_task_request_queue,
                                                       hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("select_count_distinct_jid_workers_by_poolname: %s"
                     % JTM_SQL["select_count_distinct_jid_workers_by_poolname"]
                     % dict(pool_name=inner_task_request_queue,
                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        db.close()

        pool_size = num_nodes_to_request * num_workers_per_node
        num_live_worker_in_pool = num_slurm_jid_in_pool * num_workers_per_node
        num_worker_to_add = int(ceil(float(pool_size - num_live_worker_in_pool) / num_workers_per_node))

        logger.debug("num_worker_to_add={} pool_size={} "
                     "num_live_worker_in_pool={} num_slurm_jid_in_pool={} "
                     "NUM_TOTAL_WORKERS={} num_workers_per_node={}""".format(num_worker_to_add,
                                                                             pool_size,
                                                                             num_live_worker_in_pool,
                                                                             num_slurm_jid_in_pool,
                                                                             NUM_TOTAL_WORKERS.value,
                                                                             num_workers_per_node))

        uniq_worker_id = None
        for i in range(0, num_worker_to_add):
            b_failed_to_request_worker = False
            # NOTE: if pool_name is None, JTM will spawn worker(s) in the main pool, not in a custom pool
            # sbatch_cmd_str = """ssh -t {}@{}.nersc.gov "{} && jtm-worker -wt dynamic {} -cl {} -c {} -t {} -m {}" """
            # Now only consider jtm running from nersc nodes and multiple jtm instances per each platform
            #
            # NOTE: User can request only "dynamic" workers from WDL. The "static" workers are managed
            #  by the admin.
            uniq_worker_id = str(shortuuid.uuid())

            # sbatch_cmd_str = """{}jtm-worker -wt dynamic {} -cl {} -c {} -t {} \
            # -m {} -wi {} -C {} -nw {} --qos {} --account {}\
            sbatch_cmd_str = """{}jtm worker -wt dynamic -p {} -cl {} -c {} -t {} \
            -m {} -wi {} -C {} -nw {} --qos {} -A {}\
            """.format("%s && " if ENV_ACTIVATION is not None else "",
                       pool_name,
                       pool_cluster,
                       pool_ncpus,
                       pool_time,
                       pool_mem,
                       uniq_worker_id,
                       pool_constraint,
                       num_workers_per_node,
                       pool_qos,
                       pool_charge_account)

            logger.info("Executing {}".format(sbatch_cmd_str))
            so, _, ec = run_sh_command(sbatch_cmd_str, live=True, log=logger)

            # Print job script for logging
            run_sh_command(sbatch_cmd_str + " --dry-run", live=True, log=logger)

            # Get the slurm job id returned from jtm-worker
            try:
                slurm_job_id = int(so.split('\n')[1])
            except Exception:
                logger.critical("Failed to get a valid job ID back from requesting a dynamic worker")
                ec = 1  # make it fail

            if ec == 0:  # if sbatch by jtm-worker done successfully
                logger.debug("Insert into workers table.")
                for nwpn in range(num_workers_per_node):
                    # Insert into workers for new worker id
                    try:
                        db = DbSqlMy(db=MYSQL_DB)
                        # If a dynamic worker is requested successfully,
                        # insert the info into workers table
                        # loggger.info("Try to update workers table")check_worker

                        # Fixme: After sbatch, another sbatch with the same pool name will be executed
                        #  again.
                        # Solution: Set the lifeleft=-2 and update select_count_workers_by_poolname sql
                        # statement to check only lifeleft!=-1 so that num_live_worker_in_pool can include
                        # already sbatched workers for the pool.
                        #
                        db.execute(JTM_SQL["insert_workers_workerid_workertype_poolname"]
                                   % dict(worker_id=uniq_worker_id + str(nwpn+1),
                                          worker_type=WORKER_TYPE["dynamic"],
                                          pool_name=inner_task_request_queue,
                                          jtm_host_name=pool_cluster,
                                          lifeleft=-2,
                                          slurmjobid=slurm_job_id,
                                          # nwpn=num_workers_per_node,
                                          nwpn=1))
                        db.commit()
                        db.close()
                    except Exception as e:
                        logger.critical(e)
                        logger.critical("Failed to insert workers table for workerid and workertype.")
                        logger.debug("Retry to insert workers table for workerid and workertype.")
                        raise
            else:
                logger.critical("Failed to execute the command, %s" % (sbatch_cmd_str))
                logger.critical("Failed to request workers.")
                send_msg_callback(ch, method, props, TASK_STATUS["invalidtask"],
                                  JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, delivery_mode=2)
                b_failed_to_request_worker = True
                break

    if not b_failed_to_request_worker:
        lastTid = -1

        # Fixme: mysql connection.py IndexError: bytearray index out of range
        # note: seems like mysql pool connection error
        success = False
        while success is not True:
            success = True
            try:
                db = DbSqlMy(db=MYSQL_DB)
                # table fields: userCmd, outFiles, doneFlag, retryCnt, task_type
                db.execute(JTM_SQL["insert_tasks_usercmd_outfiles"]
                           % (user_task_cmd, output_file, "0", 0, TASK_TYPE[task_type]))
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
                           % dict(task_id=lastTid,
                                  status_id=TASK_STATUS["ready"]))
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
            task_status_int = int(db.selectScalar(JTM_SQL["select_status_runs_by_taskid"]
                                                  % dict(task_id=lastTid)))
            b_already_canceled = int(db.selectScalar(JTM_SQL["select_cancelled_runs_by_taskid"]
                                                     % dict(task_id=lastTid)))
            db.close()

            # Note: this is just in case
            #  It is based on the assumption that a task can be cancelled between ready -> queued
            #  status change.
            if b_already_canceled == 1 or task_status_int == TASK_STATUS["terminated"]:
                if task_status_int != TASK_STATUS["terminated"]:
                    db = DbSqlMy(db=MYSQL_DB)
                    db.execute(JTM_SQL["update_tasks_doneflag_by_taskid"]
                               % dict(task_id=lastTid,
                                      done_flag=TASK_STATUS["terminated"]))
                    db.commit()
                    db.close()

            else:
                # Prepare msg to jtm-worker
                msg_to_send_dict = {}
                msg_to_send_dict["task_id"] = lastTid
                msg_to_send_dict["user_cmd"] = user_task_cmd
                msg_to_send_dict["output_files"] = output_file
                msg_to_send_dict["done_flag"] = 0
                msg_to_send_dict["task_type"] = TASK_TYPE[task_type]
                msg_to_send_dict["output_dir"] = output_dir
                msg_to_send_dict["stdout"] = stdout_file
                msg_to_send_dict["stderr"] = stderr_file
                # msg_to_send_dict["cromwell_jid"] = cromwellJid

                logger.info("Total number of workers (alive + requested): %d", NUM_TOTAL_WORKERS.value)

                # Create and send request message to workers
                msg_zipped = zdumps(json.dumps(msg_to_send_dict))
                corr_id = str(uuid.uuid4())

                logger.info("Send a task to {}".format(inner_task_request_queue))

                try:
                    ch.basic_publish(exchange=JTM_INNER_MAIN_EXCH,
                                     routing_key=inner_task_request_queue,
                                     properties=pika.BasicProperties(
                                         delivery_mode=2,  # make message persistent
                                         reply_to=JTM_INNER_RESULT_Q,  # set reply queue name
                                         correlation_id=corr_id),
                                     body=msg_zipped)
                except Exception as detail:
                    logger.exception("Exception: Failed to send a request to %s", inner_task_request_queue)
                    logger.exception("Detail: %s", str(detail))
                    # Todo: set the task status --> failed
                    # ch.stop_consuming()
                    # ch.close()
                    # sys.exit(1)
                    raise OSError(2, 'Failed to send a request to a worker')

                # Update status to "queued"
                # Todo: need this update to change the task status to "queued"?
                time.sleep(TASK_STAT_UPDATE_INTERVAL)

                try:
                    db = DbSqlMy(db=MYSQL_DB)
                    db.execute(JTM_SQL["update_runs_tid_startdate_by_tid"]
                               % dict(status_id=TASK_STATUS["queued"],
                                      now=time.strftime("%Y-%m-%d %H:%M:%S"),
                                      task_id=lastTid))
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
        logger.debug("Sending task id = %d to jtm-submit via %s" % (lastTid, props.reply_to))
        send_msg_callback(ch, method, props, lastTid, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, delivery_mode=2)


# -------------------------------------------------------------------------------
def process_task_status(ch, method, props, task_id):
    """
    Process task status request
    :param ch: channel
    :param method: method_frame
    :param props: pika connection property
    :param task_id: for task ID
    :return:
    """

    db = DbSqlMy(db=MYSQL_DB)
    # cur = db.execute(JTM_SQL["select_status_runs_by_taskid"] % dict(task_id=task_id))
    cur = db.execute(JTM_SQL["select_status_cancelled_runs_by_taskid"]
                     % dict(task_id=task_id))
    ret = cur.fetchone()
    db.close()

    if ret is not None:
        logger.debug("Task status from runs table: {}".format(ret))
        if ret[1] == 1:  # cancellation requested
            task_status_int = -4
        else:
            task_status_int = ret[0]
    else:
        logger.debug("Task ID not found: {}".format(task_id))
        task_status_int = -5  # invalid task id

    # Send tid to jtm_status
    logger.info("Sent reply: task status = %s via %s" % (str(task_status_int), props.reply_to))
    send_msg_callback(ch, method, props, task_status_int, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, delivery_mode=2)


# -------------------------------------------------------------------------------
def process_resource_log(ch, method, props, task_id):
    """
    With given task id, get the resource log file and send to jtm-resource-log

    :param ch:
    :param method:
    :param props:
    :param task_id:
    :return:
    """
    resource_log_file = None
    db = DbSqlMy(db=MYSQL_DB)
    cur = db.execute(JTM_SQL["select_resource_runs_by_taskid"]
                     % dict(task_id=task_id))
    ret = cur.fetchone()
    db.close()
    if ret:
        resource_log_file = ret[0]

    # Send tid to jtm_status
    logger.info("Sent reply: %s via %s" % (resource_log_file, props.reply_to))
    send_msg_callback(ch, method, props, resource_log_file, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, delivery_mode=2)


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
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()

    exch = JTM_TASK_KILL_EXCH
    queue_name = JTM_TASK_KILL_Q

    ch.exchange_declare(exchange=exch,
                        exchange_type="fanout",
                        durable=True,
                        auto_delete=False)

    message = {"task_id": task_id,
               "worker_id": wid,
               "child_pid": cpid}
    msg_zipped = zdumps(json.dumps(message))

    ch.queue_declare(queue=queue_name,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)

    # Create routing key with "wid"."task_id" so that workers can filter the messages
    routing_key = str(wid)
    ch.queue_bind(exchange=exch,
                  queue=queue_name,
                  routing_key=routing_key)

    try:
        logger.info("Send send_task_kill_request to worker %s, %r with routing key %s" %
                    (queue_name, message, routing_key))
        assert queue_name.endswith(CNAME)
        ch.basic_publish(exchange=exch,
                         routing_key=queue_name,
                         properties=pika.BasicProperties(delivery_mode=2),  # make message persistent
                         body=msg_zipped)
    except Exception as e:
        logger.critical("Something wrong in send_task_kill_request(): %s", e)
        raise

    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def send_reproduce_or_die(worker_id, child_proc_id, nClones):
    """
    Let a worker kill itself or clone itself
     - static worker: nClones = 1
     - dynamic worker: nClones >= 2

    * Note: dynamic workers are no longer cloned. Thus, this function is only for cloning static
            workers.

    :param worker_id:  worker id
    :param child_proc_id: child process id
    :param nClones: number of clones to spawn (if -1, kill the worker)
    :return:
    """
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()
    ch.exchange_declare(exchange=JTM_WORKER_POISON_EXCH,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)

    poison_queue_name = JTM_WORKER_POISON_Q
    message = {"worker_id": worker_id,
               "child_pid": child_proc_id,
               "num_clones": nClones}
    msg_zipped = zdumps(json.dumps(message))

    ch.queue_declare(queue=poison_queue_name,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)
    ch.queue_bind(exchange=JTM_WORKER_POISON_EXCH,
                  queue=poison_queue_name,
                  routing_key=worker_id)

    try:
        if nClones == -1:
            logger.info("Send poison to worker %s, %r" % (poison_queue_name, message))
        else:
            logger.info("Send cloning signal to worker %s, %r" % (poison_queue_name, message))

        assert poison_queue_name.endswith(CNAME)
        ch.basic_publish(exchange=JTM_WORKER_POISON_EXCH,
                         routing_key=worker_id,
                         properties=pika.BasicProperties(delivery_mode=2),  # make message persistent
                         body=msg_zipped)
        logger.debug("send_reproduce_or_die published {} with {}".format(message, worker_id))

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
    task_id = int(msg["task_id"])
    return_msg = None

    db = DbSqlMy(db=MYSQL_DB)
    try:
        task_status = db.selectScalar(JTM_SQL["select_status_runs_by_taskid"]
                                      % dict(task_id=task_id))
    except AssertionError:
        task_status = None
    db.close()

    def update_runs_cancelled(task_id):
        db = DbSqlMy(db=MYSQL_DB)
        db.execute(JTM_SQL["update_runs_cancelled_by_tid"]
                   % dict(task_id=task_id, now=time.strftime("%Y-%m-%d %H:%M:%S")))
        db.commit()
        db.close()
        # Note: why I can't just set the "status" field to -4 = terminated here?

    if task_status:
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
        if task_status in (TASK_STATUS["ready"], TASK_STATUS["queued"]):
            logger.debug("Task cancellation requested but the task has already been queued.")
            logger.debug("The task will be terminated once it's started.")

            # Update runs table with cancellation requested (cancelled = 1)
            update_runs_cancelled(task_id)

            # Note: if jtm_status is executed, this "cancelled" field should be checked to
            #  determine the task status

            return_msg = 0
        elif task_status == TASK_STATUS["running"]:
            logger.debug("Task cancellation requested. The task is being terminated.")
            # # Now worker id and child pid is ready
            # db = DbSqlMy(db=MYSQL_DB)
            # worker_id_list = db.selectAs1Col(JTM_SQL["select_workerid_workers_by_tid"] % dict(task_id=task_id))
            # child_proc_id = db.selectScalar(JTM_SQL["select_chilpid_runs_by_tid"] % dict(task_id=task_id))
            # logger.debug("SQL: select_chilpid_runs_by_tid = %d" % (int(child_proc_id) if child_proc_id else 0))
            # logger.debug("SQL: select_workerid_workers_by_tid %d = %s" % (task_id, worker_id_list))
            # db.close()
            #
            # assert child_proc_id > 0
            # assert worker_id_list is not None
            #
            # send_task_kill_request(task_id, worker_id_list[0], child_proc_id)

            # Update runs table with cancellation requested (cancelled = 1)
            update_runs_cancelled(task_id)
            return_msg = 0
        elif task_status == TASK_STATUS["success"]:
            logger.debug("Task cancellation requested but the task is in completed status.")
            return_msg = 0
        elif task_status == TASK_STATUS["terminated"]:
            logger.debug("Task cancellation requested but the task is already in terminated status.")
            return_msg = 0
        elif task_status == TASK_STATUS["failed"]:
            logger.debug("Task cancellation requested but the task is already in failed status.")
            return_msg = 0
        elif task_status in (TASK_STATUS["outputerror"], TASK_STATUS["outputerror"], TASK_STATUS["outputerror"]):
            logger.debug("Task cancellation request is ignored.")
            return_msg = 0
        else:
            logger.debug("Failed to cancel a task. Unexpected condition.")
            return_msg = -1
    else:  # task id not found
        return_msg = -5

    # Send status to jtm_kill
    logger.debug("Return kill command result, %d to jtm-kill." % (return_msg))
    send_msg_callback(ch, method, props, return_msg, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, delivery_mode=2)


# -------------------------------------------------------------------------------
def process_check_worker(ch, method, props, msg_unzipped):
    # If custom_pool name is set, get the number of workers in the pool
    # else the total number of live workers will be sent
    logger.debug("jtm-check-worker: %s" % str(msg_unzipped))

    if "task_pool" in msg_unzipped and msg_unzipped["task_pool"] and \
            "jtm_host_name" in msg_unzipped and msg_unzipped["jtm_host_name"]:
        db = DbSqlMy(db=MYSQL_DB)
        # Try to select "pool_name" in workers table by hostname + uselifeLeftrname + poolname
        # Check timediff(now()-end_datetime) in workers table to filter out dead worker

        new_pool_name = JTM_INNER_REQUEST_Q + '.' + msg_unzipped["task_pool"]
        num_live_worker_in_pool = db.selectScalar(JTM_SQL["select_count_workers_by_jtm_host_name_poolname_enddate"]
                                                  % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                                                         pool_name=new_pool_name,
                                                         hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("select_count_workers_by_jtm_host_name_poolname_enddate: %s"
                     % JTM_SQL["select_count_workers_by_jtm_host_name_poolname_enddate"]
                     % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                            pool_name=new_pool_name,
                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("node cnt in the pool: %s" % num_live_worker_in_pool)

        num_total_num_workers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_jtm_host_name_poolname_enddate"]
                                                % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                                                       pool_name=new_pool_name,
                                                       hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("select_sum_nwpn_workers_by_jtm_host_name_poolname_enddate: %s"
                     % JTM_SQL["select_sum_nwpn_workers_by_jtm_host_name_poolname_enddate"]
                     % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                            pool_name=new_pool_name,
                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("worker cnt in the pool: %s" % num_total_num_workers)

        db.close()

        # send_msg_callback(ch, method, props, num_live_worker_in_pool, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)
        send_msg_callback(ch, method, props,
                          num_total_num_workers if num_total_num_workers else 0,
                          JGI_JTM_MAIN_EXCH,
                          JTM_TASK_RESULT_Q)

    elif "task_pool" in msg_unzipped and msg_unzipped["task_pool"]:
        db = DbSqlMy(db=MYSQL_DB)
        # Try to select "pool_name" in workers table by hostname + uselifeLeftrname + poolname
        # Check timediff(now()-end_datetime) in workers table to filter out dead worker

        new_pool_name = JTM_INNER_REQUEST_Q + '.' + msg_unzipped["task_pool"]
        num_live_worker_in_pool = db.selectScalar(JTM_SQL["select_count_workers_by_poolname_enddate"]
                                                  % dict(pool_name=new_pool_name,
                                                         hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("select_count_workers_by_poolname_enddate: %s"
                     % JTM_SQL["select_count_workers_by_poolname_enddate"]
                     % dict(pool_name=new_pool_name,
                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("node cnt in the pool: %s" % num_live_worker_in_pool)

        num_total_num_workers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_poolname_enddate"]
                                                % dict(pool_name=new_pool_name,
                                                       hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("select_sum_nwpn_workers_by_poolname_enddate: %s"
                     % JTM_SQL["select_sum_nwpn_workers_by_poolname_enddate"]
                     % dict(pool_name=new_pool_name,
                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("worker cnt in the pool: %s" % num_total_num_workers)

        db.close()

        # send_msg_callback(ch, method, props, num_live_worker_in_pool, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)
        send_msg_callback(ch, method, props,
                          num_total_num_workers if num_total_num_workers else 0,
                          JGI_JTM_MAIN_EXCH,
                          JTM_TASK_RESULT_Q)

    elif "jtm_host_name" in msg_unzipped and msg_unzipped["jtm_host_name"]:
        db = DbSqlMy(db=MYSQL_DB)
        # Check timediff(now()-end_datetime) in workers table to filter out dead worker
        num_live_workers = db.selectScalar(JTM_SQL["select_count_workers_by_jtm_host_name"]
                                           % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                                                  hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug(JTM_SQL["select_count_workers_by_jtm_host_name"]
                     % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("node cnt in the host: %s" % num_live_workers)

        num_total_num_workers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_jtm_host_name_enddate"]
                                                % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                                                       hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug(JTM_SQL["select_sum_nwpn_workers_by_jtm_host_name_enddate"]
                     % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                            hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("worker cnt in the host: %s" % num_total_num_workers)

        db.close()

        # send_msg_callback(ch, method, props, num_live_workers, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)
        send_msg_callback(ch, method, props,
                          num_total_num_workers if num_total_num_workers else 0,
                          JGI_JTM_MAIN_EXCH,
                          JTM_TASK_RESULT_Q)
    else:
        send_msg_callback(ch, method, props,
                          NUM_TOTAL_WORKERS.value,
                          JGI_JTM_MAIN_EXCH,
                          JTM_TASK_RESULT_Q)


# -------------------------------------------------------------------------------
def process_remove_pool(ch, method, props, msg_unzipped):
    task_pool_name = JTM_INNER_REQUEST_Q + '.' + msg_unzipped["task_pool"]
    db = DbSqlMy(db=MYSQL_DB)
    logger.debug(JTM_SQL["select_jid_workers_by_poolname"]
                 % dict(pool_name=task_pool_name))

    # Get the list of slurm job id to cancel
    zombie_slurm_job_id_list = db.selectAll(JTM_SQL["select_jid_workers_by_poolname"]
                                            % dict(pool_name=task_pool_name))
    for jid in zombie_slurm_job_id_list:
        # print jid[0]
        scancel_cmd = "scancel %s" % (jid[0])
        _, _, ec = run_sh_command(scancel_cmd, live=True, log=logger)
        if ec == 0:
            logger.info("Successfully cancel the job, %s" % (jid[0]))
        else:
            logger.debug("%s not found." % (jid[0]))

    # Update workers table for canceled job with life_left=-2
    # Note: the end_datetime MUST BE updated with now() so that it is not counted as alive workers
    logger.debug(JTM_SQL["update_lifeleft_enddate_workers_by_poolname"]
                 % dict(pool_name=task_pool_name,
                        now=time.strftime("%Y-%m-%d %H:%M:%S")))
    db.execute(JTM_SQL["update_lifeleft_enddate_workers_by_poolname"]
               % dict(pool_name=task_pool_name,
                      now=time.strftime("%Y-%m-%d %H:%M:%S")))
    db.close()
    logger.info("Pool, %s removed!" % task_pool_name)
    send_msg_callback(ch, method, props, 1, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)


# -------------------------------------------------------------------------------
def on_task_request(ch, method, props, body):
    """
    Event handler for processing request from jtm-* CLI tools
    """
    # Uncompress msg
    msg_unzipped = json.loads(zloads(body))

    logger.info("New task: {}".format(msg_unzipped))
    assert "task_type" in msg_unzipped, "Critical: need a task type!"
    task_type = msg_unzipped["task_type"]
    assert task_type in TASK_TYPE, "Critical: invalid task type: {}".format(msg_unzipped)
    inner_task_request_queue = JTM_INNER_REQUEST_Q

    # If pool is set, the tasks which use the queue name (=pool name)
    # will only be sent to the pool of workers
    # NOTE: pool_name might be None even if "pool" in msg_unzipped
    if "pool" in msg_unzipped:
        pool_spec_json_str = json.loads(json.dumps(msg_unzipped["pool"]))
        pool_name = pool_spec_json_str["name"] if "name" in pool_spec_json_str and pool_spec_json_str["name"] else None
        if pool_name:
            # _jtm_inner_request_queue.<cluster_name>.jtm.<pool_name>
            inner_task_request_queue = JTM_INNER_REQUEST_Q + "." + pool_name
        else:
            inner_task_request_queue = JTM_INNER_REQUEST_Q + ".small"

    # This allows jtm to keep task requests without any alive worker
    ch.exchange_declare(exchange=JTM_INNER_MAIN_EXCH,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)

    # Get queue length
    # taskQueueLen = ch.queue_declare(queue=inner_task_request_queue,
    #                                 durable=True,
    #                                 exclusive=False,
    #                                 auto_delete=True).method.message_count
    # logger.debug("#tasks queued = %d", taskQueueLen)

    ch.queue_declare(queue=inner_task_request_queue,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)
    ch.queue_bind(exchange=JTM_INNER_MAIN_EXCH,
                  queue=inner_task_request_queue,
                  routing_key=inner_task_request_queue)

    if task_type == "task":
        # Todo: inner_task_request_queue can be customized for different types of workers ==> done
        # ch.queue_bind() needs to be called with new queue name
        process_task_request(ch, method, props, msg_unzipped, inner_task_request_queue)

    elif task_type == "status":
        process_task_status(ch, method, props, int(msg_unzipped["task_id"]))

    elif task_type == "resource":
        # Note: this just returns the resource log file name. The file content should be converted
        #  into JSON format on the user side b/c the file might only be accessible from the user's
        #  machine.
        process_resource_log(ch, method, props, int(msg_unzipped["task_id"]))

    elif task_type == "kill":
        # Todo: need to check msg_unzipped["jtm_host_name"] to determine a way to kill a task per
        #  cluster/cloud
        process_task_kill(ch, method, props, msg_unzipped)

    elif task_type == "check_manager":
        # Just reply back to jtm-check-manager with '88'
        send_msg_callback(ch, method, props, 88, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)

    elif task_type == "check_worker":
        process_check_worker(ch, method, props, msg_unzipped)

    elif task_type == "remove_pool":
        # Todo: need to check msg_unzipped["jtm_host_name"] to determine a way to remove pool per
        #  cluster/cloud
        process_remove_pool(ch, method, props, msg_unzipped)

    else:
        logger.critical("Undefined task type: {}".format(task_type))
        send_msg_callback(ch, method, props, -1, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)
        logger.debug("sent reply: %s", str(-1))


# -------------------------------------------------------------------------------
# Not used
# def kill_all_workers(pool):
#     """
#     Send termination signal to all the workers
#     :param pool: worker pool name (= taskqueue name) to kill
#     :return:
#     """
#     # Remote broker (rmq.nersc.gov) connection open
#     rmq_conn = RmqConnectionHB()
#     conn = rmq_conn.open()
#     ch = conn.channel()
#     ch.exchange_declare(exchange=JTM_INNER_MAIN_EXCH,
#                         exchange_type="fanout",
#                         durable=False,
#                         auto_delete=False)
#
#     msg_to_send_dict = {}
#     msg_to_send_dict["task_type"] = TASK_TYPE["term"]
#     msg_to_send_dict["task_queue"] = pool
#     msg_zipped = zdumps(json.dumps(msg_to_send_dict))
#
#     assert pool.endswith(CNAME)
#     try:
#         ch.basic_publish(exchange=JTM_INNER_MAIN_EXCH,
#                          routing_key=pool,
#                          body=msg_zipped)
#     except Exception as detail:
#         logger.exception("Exception: Failed to submit a TERM signal to the workers.")
#         logger.exception("Detail: %s", str(detail))
#         # ch.stop_consuming()
#         # ch.close()
#         # conn.close()
#         # sys.exit(1)
#         raise
#
#     ch.close()
#     conn.close()


# -------------------------------------------------------------------------------
def task_kill_proc():
    """
    Check "runs" table if any row with cancelled request (cancelled=1).
    Check if cancelled=1 and wid!=None and status!=-4(which means it's in running status)
    Then, send a kill request to the worker with task ID.

    """
    while 1:
        db = DbSqlMy(db=MYSQL_DB)
        # Get a list of task ids where cancelled -> requested but status != terminated
        task_id_list = [int(i) for i in db.selectAs1Col(JTM_SQL["select_tids_runs_by_cancelled_and_wid"])]
        db.close()

        for tid in task_id_list:
            try:
                # Get the wid and child pid
                db = DbSqlMy(db=MYSQL_DB)
                worker_id_list = db.selectAs1Col(JTM_SQL["select_workerid_workers_by_tid"]
                                                 % dict(task_id=tid))
                child_proc_id = db.selectScalar(JTM_SQL["select_chilpid_runs_by_tid"]
                                                % dict(task_id=tid))
                logger.debug("SQL: select_chilpid_runs_by_tid = %d" % int(child_proc_id) if child_proc_id else 0)
                logger.debug("SQL: select_workerid_workers_by_tid %d = %s" % (tid, worker_id_list))
                db.close()

                assert child_proc_id > 0
                assert worker_id_list is not None

                # Send task id and process id to worker id
                send_task_kill_request(tid, worker_id_list[0], child_proc_id)

                # Update runs table with "terminated" status
                db = DbSqlMy(db=MYSQL_DB)
                db.execute(JTM_SQL["update_runs_status_to_terminated_by_tid"]
                           % dict(task_id=tid,
                                  status_id=TASK_STATUS["terminated"],
                                  now=time.strftime("%Y-%m-%d %H:%M:%S")))

                # Note: do we have to update tasks with status = "terminated"?
                #  HB recv will check if it's terminated and update tasks table properly!
                #
                db.commit()
                db.close()
            except Exception:
                logger.exception("Something goes wrong in task_kill_proc().")
                raise

        time.sleep(TASK_KILL_INTERVAL)


# -------------------------------------------------------------------------------
def zombie_worker_cleanup_proc():
    """
    Try to find cancelled slurm job
    If any, update workers table

    Todo: Very slurm dependent. Do we need this?

    """
    while 1:
        db = DbSqlMy(db=MYSQL_DB)
        slurm_job_id = [int(i) for i in db.selectAs1Col(JTM_SQL["select_slurmjid_workers_by_lifeleft"])]
        db.close()
        if len(slurm_job_id) > 0:
            logger.info("Worker checking for %s" % str(slurm_job_id))
        for j in slurm_job_id:
            # Note: slurm dependent code!
            cmd = "sacct -j %d" % j
            so, _, ec = run_sh_command(cmd, live=True, log=logger)
            if ec == 0 and so.find("CANCELLED") != -1:
                db = DbSqlMy(db=MYSQL_DB)
                db.execute(JTM_SQL["update_workers_lifeleft_by_slurmjid"]
                           % dict(slurm_jid=j,
                                  now=time.strftime("%Y-%m-%d %H:%M:%S")))
                db.commit()
                db.close()
            time.sleep(1)

        time.sleep(WORKER_KILL_INTERVAL)


# -------------------------------------------------------------------------------
def check_processes(pid_list):
    """
    Checking if the total number of processes from the manager is NUM_MANAGER_PROCS
    if not, terminate the all the proc ids

    """
    # ps_cmd = "ps -aef | grep '%s' | grep %s | grep -v grep | awk '{print $2}'" \
    #          % ("jtm manager", getpass.getuser())
    # if sys.platform.lower() == "darwin":
    #     ps_cmd = "ps -aef | grep '%s' | grep -v grep | awk '{print $2}'" \
    #              % ("jtm manager")
    # while 1:
    #     pid_list = back_ticks(ps_cmd, shell=True)
    #     if type(pid_list) is bytes:
    #         pid_list = pid_list.decode()
    #     pid_list = [int(i) for i in pid_list.split('\n') if i]
    #     logger.debug("ps_cmd: {}".format(ps_cmd))
    #     logger.debug("pid list = {}".format(pid_list))
    #     logger.debug("len(pid) = {}".format(len(pid_list)))
    #     if len(pid_list) != NUM_MANAGER_PROCS:
    #         logger.critical("JTM manager child process error.")
    #         [os.kill(i, signal.SIGTERM) for i in pid_list]
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
def manager(custom_log_dir_name: str, b_resource_usage_log_on: bool, debug: bool) -> int:
    # Log dir setting
    log_dir_name = os.path.join(config.configparser.get("JTM", "log_dir"), "log")
    if custom_log_dir_name:
        log_dir_name = custom_log_dir_name
    make_dir(log_dir_name)

    log_level = "info"
    if debug:
        log_level = "debug"

    setup_custom_logger(log_level, log_dir_name, 1, 1)

    logger.info("\n*****************\nDebug mode is %s\n*****************" % ("ON" if debug else "OFF"))
    logger.info("Set jtm log file location to %s", log_dir_name)

    RMQ_HOST = config.configparser.get("RMQ", "host")
    RMQ_PORT = config.configparser.get("RMQ", "port")
    USER_NAME = config.configparser.get("SITE", "user_name")
    PRODUCTION = False
    if config.configparser.get("JTM", "run_mode") == "prod":
        PRODUCTION = True
    JGI_JTM_MAIN_EXCH = config.configparser.get("JTM", "jgi_jtm_main_exch")
    JTM_TASK_RESULT_Q = config.configparser.get("JTM", "jtm_task_result_q")
    JTM_TASK_REQUEST_Q = config.configparser.get("JTM", "jtm_task_request_q")
    WORKER_HB_Q_POSTFIX = config.configparser.get("JTM", "worker_hb_q_postfix")

    # Remote broker (rmq.nersc.gov) connection open
    rmq_conn = RmqConnectionHB()
    conn = rmq_conn.open()
    ch = conn.channel()
    # ch.confirm_delivery()  # 11192018 to test task is not discarded
    ch.exchange_declare(exchange=JGI_JTM_MAIN_EXCH,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)
    logger.info("JTM main exchange: %s", JGI_JTM_MAIN_EXCH)

    # Declare task sending queue (client --> worker)
    #
    # This queue can be declared from worker side. BUT
    # This should be done here to enable client to send the tasks to the
    # task queue. If it is not called here, the client should be run after
    # checking if the worker is already running so that we can ensure that
    # the task queue is already declared and safe to be used by the client.
    #
    try:
        # exclusive=False -> do not remove the queue even when the connection is closed.
        ch.queue_declare(queue=JTM_TASK_REQUEST_Q,
                         durable=True,
                         exclusive=False,
                         auto_delete=True)
        ch.queue_declare(queue=JTM_TASK_RESULT_Q,
                         durable=True,
                         exclusive=False,
                         auto_delete=True)
    except Exception as detail:
        logger.exception("Exception: The queue, %s is already in use.", JTM_TASK_RESULT_Q)
        logger.exception("Detail: %s", str(detail))
        # sys.exit(1)
        return 1

    # Queue binding for getting task request from JAWS
    ch.queue_bind(exchange=JGI_JTM_MAIN_EXCH,
                  queue=JTM_TASK_REQUEST_Q,
                  routing_key=JTM_TASK_REQUEST_Q)

    logger.info("JGI Task Manager, version: %s" % (VERSION))
    logger.info("RabbitMQ broker: %s", RMQ_HOST)
    logger.info("RabbitMQ port: %s", RMQ_PORT)
    logger.info("Default task queue name: %s", JTM_TASK_REQUEST_Q)
    logger.info("Default result queue name: %s", JTM_TASK_RESULT_Q)
    logger.info("Pika version: %s", pika.__version__)
    logger.info("Database server: %s", MYSQL_HOST)
    logger.info("Database port: %s", MYSQL_PORT)
    logger.info("Database user name: %s", MYSQL_USER)
    logger.info("Database name: %s", MYSQL_DB)
    logger.info("JTM user name: %s", USER_NAME)
    logger.info("\n*****************\nRun mode is %s\n*****************"
                % ("PROD" if PRODUCTION else "DEV"))

    #
    # MySQL: prepare task table
    #
    db = DbSqlMy(db=MYSQL_DB)
    # db.ddl(JTM_SQL["set_timezone"])
    db.ddl(JTM_SQL["create_database"] % MYSQL_DB)
    db.ddl(JTM_SQL["use_database"] % MYSQL_DB)
    db.ddl(JTM_SQL["create_table_tasks"])
    db.ddl(JTM_SQL["create_table_runs"])
    # db.createIndices(names=["task_id"], table="runs")
    db.ddl(JTM_SQL["create_table_workers"])
    db.close()

    logger.debug("Main pid = {}".format(PARENT_PROCESS_ID))

    pid_list = []
    pid_list.append(PARENT_PROCESS_ID)

    # Start heartbeat checking thread
    send_hb_to_worker_proc_hdl = mp.Process(target=send_hb_to_worker_proc)
    send_hb_to_worker_proc_hdl.daemon = True
    send_hb_to_worker_proc_hdl.start()
    pid_list.append(send_hb_to_worker_proc_hdl.pid)
    logger.info("Broadcasting heartbeats to workers...")

    def proc_clean():
        if send_hb_to_worker_proc_hdl:
            send_hb_to_worker_proc_hdl.terminate()
        if rect_hb_from_worker_proc_hdl:
            rect_hb_from_worker_proc_hdl.terminate()
        if recv_result_from_worker_proc_hdl:
            recv_result_from_worker_proc_hdl.terminate()
        if task_kill_proc_hdl:
            task_kill_proc_hdl.terminate()
        if worker_cleanup_proc_hdl:
            worker_cleanup_proc_hdl.terminate()
        if check_processes_hdl:
            check_processes_hdl.terminate()

    # Start heartbeat receiving thread
    worker_hb_queue_name = WORKER_HB_Q_POSTFIX
    try:
        rect_hb_from_worker_proc_hdl = mp.Process(target=recv_hb_from_worker_proc,
                                                  args=(worker_hb_queue_name,
                                                        log_dir_name,
                                                        b_resource_usage_log_on))
        rect_hb_from_worker_proc_hdl.daemon = True
        rect_hb_from_worker_proc_hdl.start()
    except Exception as err:
        logger.exception("Worker not found: {}".format(err))
        proc_clean()
        ch.stop_consuming()
        ch.close()
        conn.close()
        sys.exit(1)

    pid_list.append(rect_hb_from_worker_proc_hdl.pid)
    logger.info("Waiting for worker\'s heartbeats from %s", worker_hb_queue_name)

    # Start task kill thread
    task_kill_proc_hdl = mp.Process(target=task_kill_proc)
    task_kill_proc_hdl.daemon = True
    task_kill_proc_hdl.start()
    pid_list.append(task_kill_proc_hdl.pid)

    # Start worker cleanup thread
    worker_cleanup_proc_hdl = mp.Process(target=zombie_worker_cleanup_proc)
    worker_cleanup_proc_hdl.daemon = True
    worker_cleanup_proc_hdl.start()
    pid_list.append(worker_cleanup_proc_hdl.pid)

    # Start result receiving thread
    # listening to JTM_INNER_RESULT_Q to which all workers will send result messages
    try:
        recv_result_from_worker_proc_hdl = mp.Process(target=recv_result_from_workers_proc)
        recv_result_from_worker_proc_hdl.daemon = True
        recv_result_from_worker_proc_hdl.start()
    except Exception as err:
        logger.exception("recv_result_from_workers_proc: {}".format(err))
        proc_clean()
        ch.stop_consuming()
        ch.close()
        conn.close()
        sys.exit(1)

    pid_list.append(recv_result_from_worker_proc_hdl.pid)

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

    logger.info("Waiting for a task request from %s", JTM_TASK_REQUEST_Q)

    # basic_qos():
    # prefetch_size (int) – This field specifies the prefetch window size. The server will send a message in advance
    # if it is equal to or smaller in size than the available prefetch size (and also falls into other prefetch limits).
    # May be set to zero, meaning “no specific limit”, although other prefetch limits may still apply. The prefetch-size
    # is ignored if the no-ack option is set in the consumer.
    # prefetch_count (int) – Specifies a prefetch window in terms of whole messages. This field may be used in
    # combination with the prefetch-size field; a message will only be sent in advance if both prefetch windows (and
    # those at the channel and connection level) allow it. The prefetch-count is ignored if the no-ack option is set in
    # the consumer.
    # all_channels (bool) – Should the QoS apply to all channels
    #
    ch.basic_qos(prefetch_count=1)
    try:
        ch.basic_consume(queue=JTM_TASK_REQUEST_Q,
                         on_message_callback=on_task_request,
                         auto_ack=False)
    except Exception as e:
        logger.exception("basic_consume: {}".format(e))
        proc_clean()
        ch.stop_consuming()
        ch.close()
        conn.close()
        sys.exit(1)

    # NOTE: the below methods might cause error in consuming messages which are already queued
    # ch.basic_consume(lambda ch, method, properties, body: on_task_request(ch, method, properties, body, g_dbConnPool),
    #                  queue=JTM_TASK_REQUEST_Q,
    #                  auto_ack=False)
    # on_task_request_callback = functools.partial(on_task_request, args=(g_dbConnPool))
    # ch.basic_consume(on_task_request_callback,
    #                  JTM_TASK_REQUEST_Q,
    #                  auto_ack=False)

    # Keep consuming messages from task request queue
    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        proc_clean()
        ch.stop_consuming()

    # unreachable
    if ch:
        ch.close()
    if conn:
        conn.close()

    return 0
