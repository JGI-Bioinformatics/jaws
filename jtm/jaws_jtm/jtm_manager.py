#! /usr/bin/env python
# pylint: disable=C0111,C0103,R0205
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)

"""

JGI task manager


Example of task processing scenario
1. jtm-submiit sends a msg to "jtm_main_exchange" with "jtm_task_request_queue" tag.
2. jtm-manager listens to "jtm_task_request_queue" which is bound to "jtm_main_exchange"
3. When a task is published, jtm-manager takes it and sends it to a pool
   (to jgi_jtm_inner_main_exchange)
4. Workers listen to jgi_jtm_inner_request_queue which is bound to "jgi_jtm_inner_main_exchange"
5. When a task is completed, a worker sends a result msg to "jgi_jtm_inner_main_exchange" with
   "jgi_jtm_inner_result_queue" tag
6. jtm-manager listens to "jgi_jtm_inner_result_queue" queue. When a result is ready,
   takes and updates tables

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
import signal

from jaws_jtm.common import setup_custom_logger, logger
from jaws_jtm.lib.sqlstmt import JTM_SQL
from jaws_jtm.lib.rabbitmqconnection import RmqConnectionHB, send_msg_callback
from jaws_jtm.lib.dbutils import DbSqlMysql
from jaws_jtm.lib.run import pad_string_path, make_dir, run_sh_command
from jaws_jtm.lib.msgcompress import zdumps, zloads


# --------------------------------------------------------------------------------------------------
# Globals
# --------------------------------------------------------------------------------------------------
NUM_TOTAL_WORKERS = mp.Value('i', 0)
DEBUG = False


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
    rmq_conn = RmqConnectionHB(config=CONFIG)
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
        raise

    # Queue binding for the queue to receive hb from workers
    ch.queue_bind(exchange=exch_name,
                  queue=hb_queue_name,
                  routing_key=hb_queue_name)

    b_is_msg_cleared = False  # all stacked messages are processed or not
    b_is_worker_found = False  # is any alive worker
    max_worker_check_count = 0  # max number of checking workers
    interval = CLIENT_HB_RECV_INTERVAL

    while True:
        worker_ids_dict = {}
        try:
            method_frame, header_frame, body = ch.basic_get(queue=hb_queue_name, auto_ack=True)
        except Exception as detail:
            logger.exception("Exception: Failed to get a message from %s.", hb_queue_name)
            logger.exception("Detail: %s", str(detail))
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

            # Worker's hb should be more than one so take all the hb's and remove redundant ids
            for i in range(int(method_frame.message_count)):
                method_frame, header_frame, body = ch.basic_get(queue=hb_queue_name, auto_ack=True)
                msg_unzipped = json.loads(zloads(body))
                msg_unzipped = {int(k): v for k, v in msg_unzipped.items()}
                a_worker_id = msg_unzipped[HB_MSG["worker_id"]]
                worker_ids_dict[a_worker_id] = msg_unzipped

            for k, v in worker_ids_dict.items():
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

                if b_resource_log:
                    logger.resource(v)

                if pool_name:
                    queue_name = JTM_INNER_REQUEST_Q + "." + pool_name
                else:
                    queue_name = ""

                success = False

                # Todo: still need to do "while" for checking db connection?
                while success is not True:
                    success = True
                    try:
                        db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                        bExists = db.selectScalar(JTM_SQL["select_exists_workers_by_workerid"]
                                                  % dict(worker_id=a_worker_id))

                        if not bExists:
                            # Todo: remove num_worker_on_the_node insertion
                            db.execute(JTM_SQL["insert_workers_workerid_slurmjobid"]
                                       % dict(worker_id=a_worker_id,
                                              slurm_jid=slurm_job_id,
                                              worker_type=worker_type,
                                              nwpn=1))
                        else:
                            # Todo: stress test needed! Test if nworkers > 1000
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
                        success = False
                        conn.sleep(1)

                if task_id > 0 and root_proc_id != child_proc_id:  # if there is a user process started
                    datetime_str = datetime.datetime.now().strftime("%Y-%m-%d")
                    if log_dest_dir:
                        log_dir_name = os.path.join(log_dest_dir, "resource")
                    else:
                        log_dir_name = "%s/resource" % (os.getcwd())

                    padded_dir_str = pad_string_path(task_id, depth=3)  # task id 226 --> 00/00/02
                    make_dir(os.path.join(log_dir_name, padded_dir_str))
                    resource_log_fname = "%s/%s/jtm_resource_%d_%s.log" \
                                         % (log_dir_name, padded_dir_str, task_id, datetime_str)

                    with open(resource_log_fname, 'a') as rf:
                        rf.write(",".join([str(i) for i in v.values()]))
                        rf.write('\n')

                    os.chmod(resource_log_fname, 0o777)

                    try:
                        # Update tasks table with "running" status == 2 if status is still 0 or 1
                        db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)

                        # Update runs table for a task
                        # Todo: really need to store full path to the log?
                        db.execute(JTM_SQL["update_runs_status_by_taskid"]
                                   % dict(status_id=TASK_STATUS["running"],
                                          task_id=task_id,
                                          worker_id=a_worker_id,
                                          child_pid=child_proc_id,
                                          resource_log=resource_log_fname),
                                   debug=False)
                        db.commit()
                        db.close()
                    except Exception as e:
                        logger.critical(e)
                        logger.critical("Failed to update runs table for status.")
                        ch.close()
                        conn.close()
                        raise

            # Set unresponsive workers as dead
            # Collect worker_ids from hb
            alive_worker_id_list = []
            for k, v in worker_ids_dict.items():
                alive_worker_id_list.append(v[HB_MSG["worker_id"]])

            # Collect worker_id which are still set as alive from workers table
            success = False
            while success is not True:
                success = True
                try:
                    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                    live_worker_id_list = db.selectAll(JTM_SQL["select_workerid_workers_by_lifeleft_jtmhostname"]
                                                       % dict(jtm_host_name=jtm_host_name),
                                                       debug=False)
                    live_worker_id_list = [i[0] for i in live_worker_id_list]

                    # Check & update workers table
                    # Set life_left as -1 for dead workers
                    for w in live_worker_id_list:
                        if w not in alive_worker_id_list:
                            # Update lifeleft --> -1
                            db.execute(JTM_SQL["update_workers_lifeleft_by_workerid"]
                                       % dict(worker_id=w,
                                              jtm_host_name=jtm_host_name))
                    db.commit()
                    db.close()
                except Exception as e:
                    logger.critical(e)
                    logger.critical("Failed to update workers table for lifeleft.")
                    logger.debug("Retry to update workers table for lifeleft.")
                    success = False
                    conn.sleep(1)

            db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
            alive_total_num_workers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_lifeleftt_jtmhostname"]
                                                      % dict(jtm_host_name=jtm_host_name),
                                                      debug=False)
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
            # intervalIncRate <= intervalIncRate * CLIENT_HB_RECEIVE_INT_INC_RATE
            if WORKER_HB_CHECK_MAX_COUNT != 0 and \
                    max_worker_check_count > WORKER_HB_CHECK_MAX_COUNT:  # hit the max checking limit
                # Close connection and kill parent and itself
                ch.close()
                conn.close()
                raise OSError(2, 'Worker not found')

            # If there no workers alive after 5 checks, set life_left to -1 for all
            if max_worker_check_count == 3:
                try:
                    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                    db.execute(JTM_SQL["update_workers_lifeleft_for_last"])
                    db.commit()
                    db.close()
                except Exception as e:
                    logger.critical(e)
                    logger.critical("Failed to update workers table for lifeleft.")
                    ch.close()
                    conn.close()
                    raise

        # Todo: Need to start with shorter interval like (0.5-1 sec) for about 30 sec
        #  from the beginning and then use the user specified interval
        #  ==> 11.13.2018 just wait CLIENT_HB_RECV_INTERVAL
        # time.sleep(interval)

        # To fix connection lost
        # Ref) https://github.com/pika/pika/issues/1224
        try:
            conn.sleep(interval)
        except Exception as e:
            logger.critical(e)
            logger.critical("RMQ connection lost.")
            ch.close()
            conn.close()
            raise

    ch.close()
    conn.close()


# -------------------------------------------------------------------------------
def recv_result_from_workers_proc():
    """
    Receive hearbeats from all the workers running
    """
    # Remote broker (mq.nersc.gov)
    rmq_conn = RmqConnectionHB(config=CONFIG)
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
        ch.close()
        conn.close()
        raise

    # Queue binding for the queue to receive hb from workers
    ch.queue_bind(exchange=exch_name,
                  queue=inner_result_queue_name,
                  routing_key=inner_result_queue_name)

    # Todo: change to a threaded version
    ch.basic_qos(prefetch_count=NUM_RESULT_RECV_THREADS)

    try:
        ch.basic_consume(queue=inner_result_queue_name,
                         on_message_callback=recv_result_on_result)
    except Exception as detail:
        logger.exception("Exception: The queue, %s is in use.", JTM_TASK_RESULT_Q)
        logger.exception("Detail: %s", str(detail))
        ch.close()
        conn.close()
        raise

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        ch.stop_consuming()
        raise

    ch.close()
    conn.close()


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
        try:
            msg_unzipped = json.loads(zloads(body))
        except Exception as e:
            logger.critical(e)
            logger.exception("Failed to load returned result message, {}".format(body))
            raise

        assert "task_id" in msg_unzipped
        task_id = int(msg_unzipped["task_id"])
        done_flag = int(msg_unzipped["done_flag"]) if "done_flag" in msg_unzipped else -2
        ret_msg = msg_unzipped["ret_msg"] if "ret_msg" in msg_unzipped else ""
        a_worker_id = msg_unzipped["worker_id"] if "worker_id" in msg_unzipped else ""
        host_name = msg_unzipped["host_name"] if "host_name" in msg_unzipped else ""

        logger.info("Result received: {}".format(msg_unzipped))

        if ret_msg != "hb":
            logger.debug("Update tasks {} with {}".format(task_id, done_flag))
            try:
                db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
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
                elif done_flag == -5:
                    task_status_int = TASK_STATUS["failed"]
                elif done_flag == -6:
                    task_status_int = TASK_STATUS["timeout"]
                elif done_flag == -7:
                    task_status_int = TASK_STATUS["lostconnection"]
                else:
                    logger.critical("Unknown return code {}".format(done_flag))
                    raise OSError(2)

            # This seems like resolving runs table lock issue
            # issue: status is not changed to 4 after the update
            ch._connection.sleep(RESULT_RECEIVE_INTERVAL)

            # Sometimes workerId2 ==> 0
            # so wait a little bit if that happened
            a_worker_id_to_check = 0
            db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
            rows = db.selectAll(JTM_SQL["select_workerid2_workers_by_wid"]
                                % dict(worker_id=a_worker_id))
            db.close()

            try:
                a_worker_id_to_check = int(rows[0][0])
            except Exception:
                a_worker_id_to_check = 0

            # Just in case
            while a_worker_id_to_check == 0:
                db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                rows = db.selectAll(JTM_SQL["select_workerid2_workers_by_wid"]
                                    % dict(worker_id=a_worker_id))
                db.close()
                try:
                    a_worker_id_to_check = int(rows[0][0])
                except Exception:
                    a_worker_id_to_check = 0

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
                    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                    db.execute(JTM_SQL["update_runs_status_workerid2_by_taskid_2"]
                               % dict(status_id=task_status_int,
                                      wid2=a_worker_id_to_check,
                                      now=time.strftime("%Y-%m-%d %H:%M:%S"),
                                      task_id=task_id),
                               debug=DEBUG)
                    db.commit()
                    logger.debug("status {}, worker id {}, taskid {}".format(task_status_int,
                                                                             a_worker_id_to_check,
                                                                             task_id))
                    db.close()

                except Exception as e:
                    logger.critical(e)
                    logger.critical("Failed to update runs table for status and workerid2.")
                    logger.debug("Retry to update runs table for status and workerid2.")
                    success = False
                    logger.debug("update_runs_status_workerid2_by_taskid_2 sleep for %f" % RUNS_INFO_UPDATE_WAIT)
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
                logger.info("Task %s --> Failed with out-of-mem. stdout = %s, worker/host: %s/%s",
                            task_id, ret_msg, a_worker_id, host_name)
            elif done_flag == DONE_FLAGS["failed with user termination"]:  # -4
                logger.info("Task %s --> Failed by user termination. stdout = %s, worker/host: %s/%s",
                            task_id, ret_msg, a_worker_id, host_name)
            elif done_flag == DONE_FLAGS["failed with input file or command not found"]:  # -5
                logger.info("Task %s --> Failed with input file or command not found. stdout = %s, worker/host: %s/%s",
                            task_id, ret_msg, a_worker_id, host_name)
            elif done_flag == DONE_FLAGS["failed with timeout"]:  # -6
                logger.info("Task %s --> Failed with timeout. stdout = %s, worker/host: %s/%s",
                            task_id, ret_msg, a_worker_id, host_name)
            else:
                logger.warning("Cannot recognize the return code: %d" % done_flag)
                logger.info("Task %s --> worker/host: %s/%s", task_id, a_worker_id, host_name)
                raise OSError(2, 'Cannot recognize the DONE_FLAG return code')
        else:
            # If it is not a result msg, return it back to the exchange
            ch.basic_reject(delivery_tag=method.delivery_tag,
                            requeue=False)

        # After this the result message will be deleted from RabbitMQ
        # If this worker crashes while running a user command, this task will
        # be sent to other workers available
        ch.basic_ack(delivery_tag=method.delivery_tag)

    else:
        logger.critical("No data found from the result message.")
        raise OSError(2, 'No data found from the result message')


# -------------------------------------------------------------------------------
def create_sbatch_cmd_lbl(pool_name,
                          pool_cluster,
                          num_cpus_to_request,
                          pool_time,
                          mem_per_node_to_request,
                          mem_per_cpu_to_request,
                          pool_partition,
                          uniq_worker_id,
                          num_workers_per_node,
                          pool_qos,
                          charging_account,
                          job_script_dir,
                          num_nodes_to_request,
                          dry_run):
    batch_job_script_file = os.path.join(job_script_dir, "jtm_dynamic_worker_%s.job"
                                         % (uniq_worker_id))
    batch_job_script_str = ""
    batch_job_misc_params = ""
    job_name = "jtm_worker_" + pool_name

    worker_config = CONFIG.config_file
    worker_config_file = CONFIG.configparser.get("JTM", "worker_config_file")
    env_activation = CONFIG.configparser.get("JTM", "env_activation")
    partition = CONFIG.configparser.get("SLURM", "partition")
    qos = CONFIG.configparser.get("SLURM", "qos")

    if not pool_partition:
        partition = pool_partition
    if not (charging_account and pool_time):
        logger.info("charging_account, jobtime: {} {}".format(charging_account,
                                                              pool_time))
        logger.critical("Need charging account, job time for LBL")
        return None
    if not pool_qos:
        qos = pool_qos

    if worker_config_file:
        worker_config = worker_config_file

    with open(batch_job_script_file, "w") as jf:
        batch_job_script_str += "#!/bin/bash -l"
        if uniq_worker_id:
            batch_job_misc_params += " -wi %(worker_id)s_${i}" \
                                     % dict(worker_id=uniq_worker_id)

        tp_param = "-p " + pool_name
        part_param = partition
        qos_param = qos
        charge_param = charging_account
        nnode_param = 1
        if num_nodes_to_request:
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
-t %(wall_time)s %(task_queue)s \
--num_worker_per_node %(num_workers_per_node)d \
-m %(mem)s \
%(other_params)s &
    sleep 1
done
wait
""" % dict(debug="--debug" if DEBUG else "",
           wall_time=pool_time,
           job_name=job_name,
           partition_name=part_param,
           qosname=qos_param,
           charging_account=charge_param,
           num_nodes_to_request=nnode_param,
           mem_per_node_setting=mnode_param,
           worker_id=uniq_worker_id,
           job_dir=job_script_dir,
           env_activation_cmd=env_activation,
           num_workers_per_node=num_workers_per_node,
           mem=mem_per_node_to_request,
           lbl_cluster_name=pool_cluster,
           worker_type="dynamic",
           task_queue=tp_param,
           other_params=batch_job_misc_params,
           export_jtm_config_file="export JTM_CONFIG_FILE=%s" % worker_config,
           set_jtm_config_file="--config=%s" % worker_config)

        if dry_run:
            logger.info(batch_job_script_str)
        else:
            jf.writelines(batch_job_script_str)

    if not dry_run:
        os.chmod(batch_job_script_file, 0o775)

    return "sbatch --parsable %s" % (batch_job_script_file)


# -------------------------------------------------------------------------------
def create_sbatch_cmd_nersc(pool_name,
                            pool_cluster,
                            num_cpus_to_request,
                            pool_time,
                            mem_per_node_to_request,
                            mem_per_cpu_to_request,
                            uniq_worker_id,
                            contraint_param,
                            num_workers_per_node,
                            qos,
                            charging_account,
                            job_script_dir,
                            num_nodes_to_request,
                            dry_run):

    batch_job_script_file = os.path.join(job_script_dir, "jtm_dynamic_worker_%s.job"
                                         % (uniq_worker_id))
    batch_job_script_str = ""
    batch_job_misc_params = ""
    job_name = "jtm_worker_" + pool_name

    worker_config = CONFIG.config_file
    worker_config_file = CONFIG.configparser.get("JTM", "worker_config_file")
    env_activation = CONFIG.configparser.get("JTM", "env_activation")
    constraint = CONFIG.configparser.get("SLURM", "constraint")
    if contraint_param:
        constraint = contraint_param

    if not (constraint and charging_account and pool_time):
        logger.info("constraint, charging_account, jobtime: {} {} {}".format(constraint,
                                                                             charging_account,
                                                                             pool_time))
        logger.critical("Need constraint, charging account, job time for Cori")
        return None

    if contraint_param:
        constraint = contraint_param

    if worker_config_file:
        worker_config = worker_config_file

    with open(batch_job_script_file, "w") as jf:
        batch_job_script_str += "#!/bin/bash -l"
        if num_nodes_to_request:
            batch_job_script_str += """
#SBATCH -N %(num_nodes_to_request)d
#SBATCH --mem=%(mem)s""" % dict(num_nodes_to_request=num_nodes_to_request,
                                mem=mem_per_node_to_request)
            batch_job_misc_params += " -N %(num_nodes_to_request)d" \
                                     % dict(num_nodes_to_request=num_nodes_to_request,)
            if num_cpus_to_request:
                batch_job_script_str += """
#SBATCH -c %(num_cores)d""" % dict(num_cores=num_cpus_to_request)
                batch_job_misc_params += " -c %(num_cores)d" \
                                         % dict(num_cores=num_cpus_to_request)

        else:
            batch_job_script_str += """
#SBATCH -c %(num_cores)d""" % dict(num_cores=num_cpus_to_request)
            batch_job_misc_params += " -c %(num_cores)d" \
                                     % dict(num_cores=num_cpus_to_request)

            if mem_per_node_to_request:
                batch_job_script_str += """
#SBATCH --mem=%(mem)s""" % dict(mem=mem_per_node_to_request)
                batch_job_misc_params += " -m %(mem)s " \
                                         % dict(mem=mem_per_node_to_request)
            else:
                batch_job_script_str += """
#SBATCH --mem-per-cpu=%(mempercore)s""" % dict(mempercore=mem_per_cpu_to_request)
                batch_job_misc_params += " -mc %(mempercore)s" \
                                         % dict(mempercore=mem_per_cpu_to_request)

            if uniq_worker_id:
                batch_job_misc_params += " -wi %(worker_id)s_${i}" \
                                         % dict(worker_id=uniq_worker_id)

        if constraint == "haswell":
            if qos:
                batch_job_script_str += """
#SBATCH -q %(qosname)s""" % dict(qosname=qos)
                batch_job_misc_params += " -q %(qosname)s" \
                                         % dict(qosname=qos)

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
            batch_job_script_str += """
#SBATCH -C knl
#SBATCH -A %(charging_account)s
#SBATCH -q %(qosname)s""" % dict(charging_account=charging_account, qosname=qos)

            batch_job_misc_params += " -A %(charging_account)s -q %(qosname)s" \
                                     % dict(charging_account=charging_account,
                                            qosname=qos)

        elif constraint == "skylake":
            batch_job_script_str += """
#SBATCH -C skylake
#SBATCH -A %(charging_account)s
#SBATCH -q %(qosname)s""" % dict(charging_account=charging_account, qosname=qos)
            batch_job_misc_params += " -A %(charging_account)s -q %(qosname)s" \
                                     % dict(charging_account=charging_account,
                                            qosname=qos)
        excl_param = ""
        if constraint != "skylake":
            excl_param = "#SBATCH --exclusive"

        tq_param = ""
        if pool_name:
            tq_param = "-p " + pool_name

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
-cl %(cluster)s \
-wt %(worker_type)s \
-t %(wall_time)s %(task_queue)s \
--num_worker_per_node %(num_workers_per_node)d \
-C %(constraint)s \
-m %(mem)s \
%(other_params)s &
    sleep 1
done
wait
    """ % dict(worker_type="dynamic",
               debug="--debug" if DEBUG else "",
               cluster=pool_cluster,
               wall_time=pool_time,
               job_dir=job_script_dir,
               worker_id=uniq_worker_id,
               task_queue=tq_param,
               num_workers_per_node=num_workers_per_node,
               env_activation_cmd=env_activation,
               other_params=batch_job_misc_params,
               constraint=constraint,
               mem=mem_per_node_to_request,
               job_name=job_name,
               exclusive=excl_param,
               export_jtm_config_file="export JTM_CONFIG_FILE=%s" % worker_config,
               set_jtm_config_file="--config=%s" % worker_config)

        if dry_run:
            logger.info(batch_job_script_str)
        else:
            jf.writelines(batch_job_script_str)

    if not dry_run:
        os.chmod(batch_job_script_file, 0o775)

    return "sbatch --parsable %s" % (batch_job_script_file)


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
    # $ jtm submit -cmd 'ls' -cl cori -p test2 -t "00:10:00" -c 1 -s 1 -m 5G

    b_failed_to_request_worker = False
    if "pool" in msg and "name" in msg["pool"] and "time" in msg["pool"]:
        pool_spec_json_str = json.loads(json.dumps(msg["pool"]))

        # Set default values defined in Config.py
        assert "name" in pool_spec_json_str
        pool_name = pool_spec_json_str["name"]
        assert pool_name is not None and pool_name != ""
        pool_cluster = POOL_CLUSTER
        pool_ncpus = POOL_NCPUS
        pool_mem = POOL_MEM
        pool_constraint = POOL_CONSTRAINT
        pool_charge_account = POOL_CHARGE_ACCOUNT
        pool_qos = POOL_QOS

        # Note: pool size = num_nodes_to_request * num_workers_per_node
        num_workers_per_node = NUM_WORKERS_PER_NODE
        num_nodes_to_request = NUM_NODES_TO_REQUEST

        # Worker type is restricted to "dynamic" for now.
        # Todo: add the feature to remove the custom pool by user task to support "static"
        #   workers in custom pool
        # NOTE: if pool_size is set, which means user sets the number of workers needed,
        # n number of jtm dynamic workers will be created, and the pool of n dynamic
        # workers won't be increased over n. The workers will be terminated if there is
        # no tasks requested
        #
        # if not,
        # only one dynamic worker will be created. The worker will be terminated if there is
        # no tasks in the queue for a specified time duration
        if "cluster" in pool_spec_json_str:  # cluster/clouds name
            pool_cluster = pool_spec_json_str["cluster"]
        if "time" in pool_spec_json_str:  # wallclocktime request
            pool_time = pool_spec_json_str["time"]
        if "cpu" in pool_spec_json_str:  # number of cores request
            pool_ncpus = int(pool_spec_json_str["cpu"])
        if "mem" in pool_spec_json_str:  # node memory request
            pool_mem = pool_spec_json_str["mem"]
        if "mempercpu" in pool_spec_json_str:  # node memory request
            pool_mempercpu = pool_spec_json_str["mempercpu"]
        if "constraint" in pool_spec_json_str:  # [haswell | knl | skylake]
            pool_constraint = pool_spec_json_str["constraint"]
        if "qos" in pool_spec_json_str:
            pool_qos = pool_spec_json_str["qos"]  # ["genepool_special", "genepool_shared", "jgi_shared", "jgi_exvivo"]
        if "account" in pool_spec_json_str:
            pool_charge_account = pool_spec_json_str["account"]  # for example, gtrqc for skylake, fungalp for the rest
        if "nwpn" in pool_spec_json_str:
            num_workers_per_node = int(pool_spec_json_str["nwpn"])
        if "node" in pool_spec_json_str:
            num_nodes_to_request = int(pool_spec_json_str["node"])
        if "partition" in pool_spec_json_str:
            pool_partition = pool_spec_json_str["partition"]
        if "dryrun" in pool_spec_json_str:
            dry_run = pool_spec_json_str["dryrun"]

        assert (len(user_task_cmd) <= 1024)
        assert (len(output_file) <= 1024)

        ################
        # Create a pool
        ################
        # Steps
        # 1. check the number if workers in the custom pool
        # 2. if 0, create a pool
        #    else send tasks to the pool

        # Todo: maintain the number of workers sbatched -->
        #   num_node_to_add = pool_size - num_live_worker_in_pool - nWorkerSbatched
        db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
        num_slurm_jid_in_pool = db.selectScalar(JTM_SQL["select_count_distinct_jid_workers_by_poolname"]
                                                % dict(pool_name=inner_task_request_queue,
                                                       hbinterval=WORKER_HB_RECV_INTERVAL * 3),
                                                debug=False)

        # TODO: slurm jid --> sacct --> doublecheck the status of node allocation and worker
        #  if no real alive workers --> request new node
        #
        #  if num_node_to_add == 0, squeue -j slurm jid
        #
        ####################################################################################
        slurm_jid_list = db.selectAll(JTM_SQL["select_distinct_jid_workers_by_poolname"]
                                      % dict(pool_name=inner_task_request_queue,
                                             hbinterval=WORKER_HB_RECV_INTERVAL * 3))
        logger.debug("slurm job id for {}: {}".format(inner_task_request_queue, slurm_jid_list))
        db.close()

        for jid in slurm_jid_list:
            cmd = "squeue -j %d" % jid
            so, _, ec = run_sh_command(cmd, log=logger, show_stdout=False)
            if ec != 0:
                num_slurm_jid_in_pool -= 1

                # Delete worker info
                logger.info("Found dead worker(s) from %s, %s" % (str(jid), so))
                db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                db.execute(JTM_SQL["delete_from_workers_by_slurmjid"]
                           % dict(slurm_jid=jid, ),
                           debug=DEBUG)
                db.commit()
                db.close()
        ####################################################################################

        pool_size = num_nodes_to_request * num_workers_per_node
        # This is the actual number of workers = jid * nwpn
        num_live_worker_in_pool = num_slurm_jid_in_pool * num_workers_per_node
        num_node_to_add = int(ceil(float(pool_size - num_live_worker_in_pool) / num_workers_per_node))

        logger.debug("num_node_to_add={} pool_size={} "
                     "num_live_worker_in_pool={} num_slurm_jid_in_pool={} "
                     "NUM_TOTAL_WORKERS={} num_workers_per_node={}""".format(num_node_to_add,
                                                                             pool_size,
                                                                             num_live_worker_in_pool,
                                                                             num_slurm_jid_in_pool,
                                                                             NUM_TOTAL_WORKERS.value,
                                                                             num_workers_per_node))

        uniq_worker_id = None
        # Job dir setting
        job_script_dir_name = os.path.join(CONFIG.configparser.get("JTM", "log_dir"), "job")
        make_dir(job_script_dir_name)

        for i in range(0, num_node_to_add):
            uniq_worker_id = str(shortuuid.uuid())
            if pool_cluster == "cori":
                sbatch_cmd_str = create_sbatch_cmd_nersc(pool_name,
                                                         pool_cluster,
                                                         pool_ncpus,
                                                         pool_time,
                                                         pool_mem,
                                                         pool_mempercpu,
                                                         uniq_worker_id,
                                                         pool_constraint,
                                                         num_workers_per_node,
                                                         pool_qos,
                                                         pool_charge_account,
                                                         job_script_dir_name,
                                                         1,
                                                         dry_run)
            elif pool_cluster == "lbl":
                sbatch_cmd_str = create_sbatch_cmd_lbl(pool_name,
                                                       pool_cluster,
                                                       pool_ncpus,
                                                       pool_time,
                                                       pool_mem,
                                                       pool_mempercpu,
                                                       pool_partition,
                                                       uniq_worker_id,
                                                       num_workers_per_node,
                                                       pool_qos,
                                                       pool_charge_account,
                                                       job_script_dir_name,
                                                       1,
                                                       dry_run)

            b_failed_to_request_worker = True
            if not dry_run and sbatch_cmd_str is not None and sbatch_cmd_str != "":
                # User can request only "dynamic" workers from WDL.
                # sbatch_cmd_str = """{}jtm {} worker \
                #     -wt dynamic \
                #     -p {} \
                #     -cl {} \
                #     -c {} \
                #     -t {} \
                #     -m {} \
                #     -wi {} {} \
                #     -nwpn {} \
                #     --qos {} \
                #     -A {}""".format("%s && " % ENV_ACTIVATION if ENV_ACTIVATION else "",
                #                     "--config=%s" % CONFIG.config_file if CONFIG else "",
                #                     pool_name,
                #                     pool_cluster,
                #                     pool_ncpus,
                #                     pool_time,
                #                     pool_mem,
                #                     uniq_worker_id,
                #                     "-C %s" % pool_constraint if pool_constraint else "",
                #                     num_workers_per_node,
                #                     pool_qos,
                #                     pool_charge_account)

                if dry_run:
                    logger.info(sbatch_cmd_str)
                    b_failed_to_request_worker = False
                else:
                    logger.info("Executing {}".format(sbatch_cmd_str))
                    # Run sbatch from jtm worker command
                    so, _, ec = run_sh_command(sbatch_cmd_str, log=logger)
                    assert ec == 0, "Failed to run sbatch commands: %s" % sbatch_cmd_str

                    # # Print job script for logging
                    # run_sh_command(sbatch_cmd_str + " --dry_run", log=logger)

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
                                db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
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
                                                  slurm_jid=slurm_job_id,
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

    if dry_run:
        send_msg_callback(ch, method, props, -87,
                          JGI_JTM_MAIN_EXCH,
                          JTM_TASK_RESULT_Q,
                          delivery_mode=2)

    elif not b_failed_to_request_worker:
        lastTid = -1

        # Fixme: mysql connection.py IndexError: bytearray index out of range
        # note: seems like mysql pool connection error
        success = False
        while success is not True:
            success = True
            try:
                db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                db.execute(JTM_SQL["insert_tasks_usercmd_outfiles"]
                           % (user_task_cmd, output_file, "0", 0, TASK_TYPE[task_type]))
                db.commit()
                lastTid = db.selectScalar(JTM_SQL["select_last_insert_id"])
                logger.debug("lastTid = %d" % lastTid)
                db.close()
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed to insert tasks table for new user task.")
                logger.debug("Retry to insert tasks table for a user task.")
                success = False
                ch._connection.sleep(1.0)

        success = False
        while success is not True:
            success = True
            try:
                db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                db.execute(JTM_SQL["insert_runs_tid_sid"]
                           % dict(task_id=lastTid,
                                  status_id=TASK_STATUS["ready"]))
                db.commit()
                db.close()
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed to insert runs table for a new run.")
                logger.debug("Retry to insert runs table for a new run.")
                success = False
                ch._connection.sleep(1.0)

        # Todo: 01282019 now every cluster has a jtm manager running, so it doesn't need to spawn
        #  dynamic worker by a static worker
        if lastTid != -1:  # if it successfully updates runs table and gets a task id
            # Todo: Check if cancelled or terminated
            db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
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
                    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
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
                    raise OSError(2, 'Failed to send a request to a worker')

                # Update status to "queued"
                # Todo: need this update to change the task status to "queued"?
                ch._connection.sleep(TASK_STAT_UPDATE_INTERVAL)

                try:
                    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                    db.execute(JTM_SQL["update_runs_tid_startdate_by_tid"]
                               % dict(status_id=TASK_STATUS["queued"],
                                      now=time.strftime("%Y-%m-%d %H:%M:%S"),
                                      task_id=lastTid),
                               debug=DEBUG)
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

    else:
        logger.debug("Sending task id = %d to jtm-submit via %s" % (-5, props.reply_to))
        send_msg_callback(ch, method, props, -5, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q, delivery_mode=2)


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

    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
    cur = db.execute(JTM_SQL["select_status_cancelled_runs_by_taskid"]
                     % dict(task_id=task_id))
    ret = cur.fetchone()
    db.close()

    if ret is not None and ret != "":
        logger.debug("Task status from runs table: status, cancelled = {}".format(ret))
        if ret[1] == 1:  # cancellation requested
            task_status_int = TASK_STATUS["terminated"]
        else:
            task_status_int = ret[0]
    else:
        logger.debug("Task ID not found: {}".format(task_id))
        task_status_int = TASK_STATUS["invalidtask"]  # invalid task id

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
    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
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
    rmq_conn = RmqConnectionHB(config=CONFIG)
    conn = rmq_conn.open()
    ch = conn.channel()

    exch = JTM_TASK_KILL_EXCH
    queue_name = JTM_TASK_KILL_Q

    # Here worker id with task id to cancel are sent to all workers
    # using "fanout".
    # Each worker will check if the worker id matches.
    # if matched, the worker kills the pid
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
        logger.info("Send send_task_kill_request to worker %s, %r with routing key %s, for task id %d" %
                    (queue_name, message, routing_key, task_id))
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
def process_task_kill(ch, method, props, msg):
    """
    Process request to terminate a task
    Steps
    # Parse the msg from jtm-kill
    # Get workerid
    # Delete tasks row for task_id (runs table will be cascaded)
    # Send kill msg

    :param ch:
    :param method:
    :param props:
    :param msg:
    :return:
    """

    task_id = int(msg["task_id"])
    return_msg = None

    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
    task_status = db.selectScalar(JTM_SQL["select_status_runs_by_taskid"]
                                  % dict(task_id=task_id))
    db.close()

    def update_runs_cancelled(task_id):
        db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
        logger.debug(db.selectAll(
            "select taskId, status, cancelled from runs where taskId = %(task_id)s"
            % dict(task_id=task_id, )))
        db.execute(JTM_SQL["update_runs_cancelled_by_tid"]
                   % dict(task_id=task_id,
                          now=time.strftime("%Y-%m-%d %H:%M:%S")),
                   debug=DEBUG)
        db.commit()
        logger.debug(db.selectAll(
            "select taskId, status, cancelled from runs where taskId = %(task_id)s"
            % dict(task_id=task_id, )))
        logger.debug("taskId, workerId2, childPid, cancelled = {}".format(
            db.selectAll("select taskId, workerId2, childPid, cancelled "
                         "from runs "
                         "where cancelled = 1 AND workerId2 != 0 AND childPid != 0")))
        db.close()
        # Note: why can't just set the "status" field to -4 = terminated here?
        # because worker id and pid is not known yet

    if task_status:
        if task_status in (TASK_STATUS["ready"], TASK_STATUS["queued"]):
            logger.debug("Task cancellation requested but the task, %d has already been queued." % (task_id))
            logger.debug("The task will be terminated once it's started.")
            update_runs_cancelled(task_id)
            return_msg = 0
        elif task_status == TASK_STATUS["running"]:
            logger.debug("Task cancellation requested. The task, %s is being terminated." % (task_id))
            update_runs_cancelled(task_id)

            # Get the wid and child pid
            db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
            worker_id_list = db.selectAs1Col(JTM_SQL["select_workerid_workers_by_tid"]
                                             % dict(task_id=task_id),
                                             debug=False)
            logger.debug("worker list = {}".format(worker_id_list))
            child_proc_id = db.selectScalar(JTM_SQL["select_chilpid_runs_by_tid"]
                                            % dict(task_id=task_id),
                                            debug=False)
            logger.debug("child_proc_id = {}".format(child_proc_id))
            db.close()

            if child_proc_id > 0 and worker_id_list is not None:
                # Send task id and process id to worker id
                logger.info("Send task kill command: {} {} {}".format(task_id, worker_id_list, child_proc_id))
                for wid in worker_id_list:
                    send_task_kill_request(task_id, wid, child_proc_id)

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
        elif task_status in (TASK_STATUS["outputerror"],
                             TASK_STATUS["outputerror"],
                             TASK_STATUS["outputerror"]):
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
        db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
        # Try to select "pool_name" in workers table by hostname + uselifeLeftrname + poolname
        # Check timediff(now()-end_datetime) in workers table to filter out dead worker

        new_pool_name = JTM_INNER_REQUEST_Q + '.' + msg_unzipped["task_pool"]
        num_live_worker_in_pool = db.selectScalar(JTM_SQL["select_count_workers_by_jtm_host_name_poolname_enddate"]
                                                  % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                                                         pool_name=new_pool_name,
                                                         hbinterval=WORKER_HB_RECV_INTERVAL * 3),
                                                  debug=False)
        logger.debug("node cnt in the pool: %s" % num_live_worker_in_pool)

        num_total_num_workers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_jtm_host_name_poolname_enddate"]
                                                % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                                                       pool_name=new_pool_name,
                                                       hbinterval=WORKER_HB_RECV_INTERVAL * 3),
                                                debug=False)
        logger.debug("worker cnt in the pool: %s" % num_total_num_workers)

        db.close()

        send_msg_callback(ch, method, props,
                          num_total_num_workers if num_total_num_workers else 0,
                          JGI_JTM_MAIN_EXCH,
                          JTM_TASK_RESULT_Q)

    elif "task_pool" in msg_unzipped and msg_unzipped["task_pool"]:
        db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
        # Try to select "pool_name" in workers table by hostname + uselifeLeftrname + poolname
        # Check timediff(now()-end_datetime) in workers table to filter out dead worker

        new_pool_name = JTM_INNER_REQUEST_Q + '.' + msg_unzipped["task_pool"]
        num_live_worker_in_pool = db.selectScalar(JTM_SQL["select_count_workers_by_poolname_enddate"]
                                                  % dict(pool_name=new_pool_name,
                                                         hbinterval=WORKER_HB_RECV_INTERVAL * 3),
                                                  debug=False)
        logger.debug("node cnt in the pool: %s" % num_live_worker_in_pool)

        num_total_num_workers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_poolname_enddate"]
                                                % dict(pool_name=new_pool_name,
                                                       hbinterval=WORKER_HB_RECV_INTERVAL * 3),
                                                debug=False)
        logger.debug("worker cnt in the pool: %s" % num_total_num_workers)
        db.close()

        send_msg_callback(ch, method, props,
                          num_total_num_workers if num_total_num_workers else 0,
                          JGI_JTM_MAIN_EXCH,
                          JTM_TASK_RESULT_Q)

    elif "jtm_host_name" in msg_unzipped and msg_unzipped["jtm_host_name"]:
        db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
        # Check timediff(now() - end_datetime) in workers table to filter out dead worker
        num_live_workers = db.selectScalar(JTM_SQL["select_count_workers_by_jtm_host_name"]
                                           % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                                                  hbinterval=WORKER_HB_RECV_INTERVAL * 3),
                                           debug=False)
        logger.debug("node cnt in the host: %s" % num_live_workers)

        num_total_num_workers = db.selectScalar(JTM_SQL["select_sum_nwpn_workers_by_jtm_host_name_enddate"]
                                                % dict(jtm_host_name=msg_unzipped["jtm_host_name"],
                                                       hbinterval=WORKER_HB_RECV_INTERVAL * 3),
                                                debug=False)
        logger.debug("worker cnt in the host: %s" % num_total_num_workers)
        db.close()

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
    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
    # Get all job id
    zombie_slurm_job_id_list = db.selectAll(JTM_SQL["select_all_jid_workers_by_poolname"]
                                            % dict(pool_name=task_pool_name,),
                                            debug=False)
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

    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
    logger.debug("process_remove_pool ****************************************** ")
    db.execute(JTM_SQL["update_runs_status_cancelled_by_status_workerId2_poolname"]
               % dict(pool_name=task_pool_name, ),
               debug=DEBUG)
    # Delete all the workers with pool name from workers table
    db.execute(JTM_SQL["delete_from_workers_by_poolname"]
               % dict(pool_name=task_pool_name,),
               debug=DEBUG)
    db.commit()
    db.close()

    logger.info("Pool, %s removed!" % task_pool_name)

    send_msg_callback(ch, method, props, 1, JGI_JTM_MAIN_EXCH, JTM_TASK_RESULT_Q)


# -------------------------------------------------------------------------------
def on_task_request(ch, method, props, body):
    """
    Event handler for processing request from jtm-* CLI tools
    """
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
            # Naming: _jtm_inner_request_queue.<cluster_name>.jtm.<pool_name>
            inner_task_request_queue = JTM_INNER_REQUEST_Q + "." + pool_name
        else:
            inner_task_request_queue = JTM_INNER_REQUEST_Q + ".small"

    # This allows jtm to keep task requests without any alive worker
    ch.exchange_declare(exchange=JTM_INNER_MAIN_EXCH,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)
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
def task_kill_proc():
    """
    Check "runs" table if any row with canceled request (canceled=1).
    Check if canceled=1 and wid!=None and status!=-4(which means it's in running status)
    Then, send a kill request to the worker with task ID.

    """
    while True:
        db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
        # Get a list of task ids where canceled -> requested but status != terminated
        tids = [int(i) for i in db.selectAs1Col(JTM_SQL["select_tids_runs_by_cancelled_and_wid2"])]
        logger.debug("tid to kill: {}".format(tids))
        logger.debug("task id list to kill: {}".format(tids))
        db.close()

        for tid in tids:
            try:
                # Get the wid and child pid
                db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                worker_id_list = db.selectAs1Col(JTM_SQL["select_workerid_workers_by_tid"]
                                                 % dict(task_id=tid),
                                                 debug=False)
                logger.debug("worker list = {}".format(worker_id_list))
                child_proc_id = db.selectScalar(JTM_SQL["select_chilpid_runs_by_tid"]
                                                % dict(task_id=tid),
                                                debug=False)
                logger.debug("child_proc_id = {}".format(child_proc_id))
                db.close()

                assert child_proc_id > 0
                assert worker_id_list is not None

                # Send task id and process id to worker id
                logger.info("Send task kill command: {} {} {}".format(tid, worker_id_list, child_proc_id))
                for wid in worker_id_list:
                    send_task_kill_request(tid, wid, child_proc_id)

                # Update runs table with "terminated" status
                logger.debug("Update runs table with terminated for tid {}".format(tid))
                db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                db.execute(JTM_SQL["update_runs_status_cancelled_to_terminated_by_tid"]
                           % dict(task_id=tid,
                                  status_id=TASK_STATUS["terminated"],
                                  now=time.strftime("%Y-%m-%d %H:%M:%S")),
                           debug=DEBUG)
                db.commit()
                db.close()
            except IndexError as e:
                logger.exception("Failed to call send_task_kill_request(): {}".format(e))
                logger.debug("tid, worker_id_list, child_proc_id: {} {} {}".format(tid,
                                                                                   worker_id_list,
                                                                                   child_proc_id))
            except Exception as e:
                logger.exception("Something goes wrong in task_kill_proc(): {}".format(e))
                raise

        time.sleep(TASK_KILL_INTERVAL)


# -------------------------------------------------------------------------------
def zombie_worker_cleanup_proc():
    """
    Try to find invalid slurm job
    If any, update workers table

    Todo: Very slurm dependent. Do we need this?

    """
    while True:
        db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
        slurm_job_id = [int(i) for i in db.selectAs1Col(JTM_SQL["select_slurmjid_workers_by_lifeleft"])]
        db.close()
        if len(slurm_job_id) > 0:
            logger.info("Worker checking for %s" % str(slurm_job_id))
        for j in slurm_job_id:
            # Note: slurm dependent code!
            cmd = "squeue -j %d" % j
            so, _, ec = run_sh_command(cmd, log=logger, show_stdout=False)
            if ec != 0:
                logger.info("Found dead worker(s) from %s, %s" % (str(j), so))
                db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
                # Note: update status from runs where workerId2 = (select workerId2 from workers)
                logger.debug("zombie_worker_cleanup_proc ****************************** ")
                db.execute(JTM_SQL["update_runs_status_cancelled_by_status_workerId2_jid"]
                           % dict(slurm_jid=j,),
                           debug=DEBUG)
                # Delete workers
                db.execute(JTM_SQL["delete_from_workers_by_slurmjid"]
                           % dict(slurm_jid=j,),
                           debug=DEBUG)
                db.commit()
                db.close()
            else:
                logger.debug("Job id {} is alive.".format(j))

            time.sleep(1)

        time.sleep(WORKER_KILL_INTERVAL)


# -------------------------------------------------------------------------------
def check_processes(pid_list):
    """
    Checking if the total number of processes from the manager is NUM_MANAGER_PROCS
    if not, terminate the all the proc ids

    """
    while True:
        logger.debug("Total Number of processes of the manager = %d" % (len(pid_list)+2))
        if len(pid_list) != NUM_MANAGER_PROCS - 2:
            raise OSError(2, 'Number of processes is wrong')
        time.sleep(NUM_PROCS_CHECK_INTERVAL)


# -------------------------------------------------------------------------------
def proc_clean(plist):
    """

    :param plist: process handle list
    :return:
    """
    for p in plist:
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
def manager(ctx: object, custom_log_dir_name: str,
            b_resource_usage_log_on: bool) -> int:

    global CONFIG
    CONFIG = ctx.obj['config']
    debug = ctx.obj['debug']
    # config file has precedence
    config_debug = CONFIG.configparser.getboolean("SITE", "debug")
    if config_debug:
        debug = config_debug
    global DEBUG
    DEBUG = debug
    global TASK_STATUS
    TASK_STATUS = CONFIG.constants.TASK_STATUS
    VERSION = CONFIG.constants.VERSION
    global TASK_TYPE
    TASK_TYPE = CONFIG.constants.TASK_TYPE
    global WORKER_TYPE
    WORKER_TYPE = CONFIG.constants.WORKER_TYPE
    global DONE_FLAGS
    DONE_FLAGS = CONFIG.constants.DONE_FLAGS
    global HB_MSG
    HB_MSG = CONFIG.constants.HB_MSG
    global NUM_MANAGER_PROCS
    NUM_MANAGER_PROCS = CONFIG.constants.NUM_MANAGER_PROCS
    global PARENT_PROCESS_ID
    PARENT_PROCESS_ID = os.getpid()  # parent process id

    MYSQL_HOST = CONFIG.configparser.get("MYSQL", "host")
    MYSQL_USER = CONFIG.configparser.get("MYSQL", "user")
    MYSQL_PORT = CONFIG.configparser.getint("MYSQL", "port")
    global MYSQL_DB
    MYSQL_DB = CONFIG.configparser.get("MYSQL", "db")
    global JTM_INNER_MAIN_EXCH
    JTM_INNER_MAIN_EXCH = CONFIG.configparser.get("JTM", "jtm_inner_main_exch")
    global CNAME
    CNAME = CONFIG.configparser.get("SITE", "instance_name")
    global JTM_TASK_RESULT_Q
    JTM_TASK_RESULT_Q = CONFIG.configparser.get("JTM", "jtm_task_result_q")
    global JTM_TASK_REQUEST_Q
    JTM_TASK_REQUEST_Q = CONFIG.configparser.get("JTM", "jtm_task_request_q")
    global JGI_JTM_MAIN_EXCH
    JGI_JTM_MAIN_EXCH = CONFIG.configparser.get("JTM", "jgi_jtm_main_exch")
    global JTM_INNER_REQUEST_Q
    JTM_INNER_REQUEST_Q = CONFIG.configparser.get("JTM", "jtm_inner_request_q")
    global JTM_INNER_RESULT_Q
    JTM_INNER_RESULT_Q = CONFIG.configparser.get("JTM", "jtm_inner_result_q")
    global WORKER_HB_RECV_INTERVAL
    WORKER_HB_RECV_INTERVAL = CONFIG.configparser.getfloat("JTM", "worker_hb_recv_interval")
    global JTM_WORKER_POISON_EXCH
    JTM_WORKER_POISON_EXCH = CONFIG.configparser.get("JTM", "jtm_worker_poison_exch")
    global JTM_WORKER_POISON_Q
    JTM_WORKER_POISON_Q = CONFIG.configparser.get("JTM", "jtm_worker_poison_q")
    global JTM_TASK_KILL_EXCH
    JTM_TASK_KILL_EXCH = CONFIG.configparser.get("JTM", "jtm_task_kill_exch")
    global JTM_TASK_KILL_Q
    JTM_TASK_KILL_Q = CONFIG.configparser.get("JTM", "jtm_task_kill_q")
    global TASK_STAT_UPDATE_INTERVAL
    TASK_STAT_UPDATE_INTERVAL = CONFIG.configparser.getfloat("JTM", "task_stat_update_interval")
    global JTM_CLIENT_HB_EXCH
    JTM_CLIENT_HB_EXCH = CONFIG.configparser.get("JTM", "jtm_client_hb_exch")
    global JTM_WORKER_HB_EXCH
    JTM_WORKER_HB_EXCH = CONFIG.configparser.get("JTM", "jtm_worker_hb_exch")
    global WORKER_KILL_INTERVAL
    WORKER_KILL_INTERVAL = CONFIG.configparser.getfloat("JTM", "worker_kill_interval")
    global TASK_KILL_INTERVAL
    TASK_KILL_INTERVAL = CONFIG.configparser.getfloat("JTM", "task_kill_interval")
    global CLIENT_HB_SEND_INTERVAL
    CLIENT_HB_SEND_INTERVAL = CONFIG.configparser.getfloat("JTM", "client_hb_send_interval")
    global NUM_RESULT_RECV_THREADS
    NUM_RESULT_RECV_THREADS = CONFIG.configparser.getint("JTM", "num_result_recv_threads")
    global NUM_PROCS_CHECK_INTERVAL
    NUM_PROCS_CHECK_INTERVAL = CONFIG.configparser.getfloat("JTM", "num_procs_check_interval")
    global ENV_ACTIVATION
    ENV_ACTIVATION = CONFIG.configparser.get("JTM", "env_activation")
    global RESULT_RECEIVE_INTERVAL
    RESULT_RECEIVE_INTERVAL = CONFIG.configparser.getfloat("JTM", "result_receive_interval")
    global RUNS_INFO_UPDATE_WAIT
    RUNS_INFO_UPDATE_WAIT = CONFIG.configparser.getfloat("JTM", "runs_info_update_wait")
    global WORKER_INFO_UPDATE_WAIT
    WORKER_INFO_UPDATE_WAIT = CONFIG.configparser.getfloat("JTM", "worker_info_update_wait")
    global CLIENT_HB_RECV_INTERVAL
    CLIENT_HB_RECV_INTERVAL = CONFIG.configparser.getfloat("JTM", "client_hb_recv_interval")
    global WORKER_HB_CHECK_MAX_COUNT
    WORKER_HB_CHECK_MAX_COUNT = CONFIG.configparser.getint("JTM", "worker_hb_check_max_count")

    # Log dir setting
    log_dir_name = os.path.join(CONFIG.configparser.get("JTM", "log_dir"), "log")
    if custom_log_dir_name:
        log_dir_name = custom_log_dir_name
    make_dir(log_dir_name)

    log_level = "info"
    if DEBUG:
        log_level = "debug"

    print("JTM Manager, version: {}".format(VERSION))

    setup_custom_logger(log_level, log_dir_name, 1, 1)

    logger.info("\n*****************\nDebug mode is %s\n*****************"
                % ("ON" if DEBUG else "OFF"))
    logger.info("Set jtm log file location to %s", log_dir_name)

    RMQ_HOST = CONFIG.configparser.get("RMQ", "host")
    RMQ_PORT = CONFIG.configparser.get("RMQ", "port")
    USER_NAME = CONFIG.configparser.get("SITE", "user_name")
    PRODUCTION = False
    if CONFIG.configparser.get("JTM", "run_mode") == "prod":
        PRODUCTION = True
    JGI_JTM_MAIN_EXCH = CONFIG.configparser.get("JTM", "jgi_jtm_main_exch")
    JTM_TASK_RESULT_Q = CONFIG.configparser.get("JTM", "jtm_task_result_q")
    JTM_TASK_REQUEST_Q = CONFIG.configparser.get("JTM", "jtm_task_request_q")
    WORKER_HB_Q_POSTFIX = CONFIG.configparser.get("JTM", "worker_hb_q_postfix")

    global POOL_CLUSTER
    POOL_CLUSTER = CONFIG.configparser.get("JTM", "cluster")
    global POOL_NCPUS
    POOL_NCPUS = CONFIG.configparser.getint("SLURM", "ncpus")
    global POOL_MEM
    POOL_MEM = CONFIG.configparser.get("SLURM", "mempernode")
    global POOL_CONSTRAINT
    POOL_CONSTRAINT = CONFIG.configparser.get("SLURM", "constraint")
    global POOL_CHARGE_ACCOUNT
    POOL_CHARGE_ACCOUNT = CONFIG.configparser.get("SLURM", "charge_accnt")
    global POOL_QOS
    POOL_QOS = CONFIG.configparser.get("SLURM", "qos")
    global NUM_WORKERS_PER_NODE
    NUM_WORKERS_PER_NODE = CONFIG.configparser.getint("JTM", "num_workers_per_node")
    global NUM_NODES_TO_REQUEST
    NUM_NODES_TO_REQUEST = CONFIG.configparser.getint("SLURM", "nnodes")

    # Remote broker (rmq.nersc.gov) connection open
    rmq_conn = RmqConnectionHB(config=CONFIG)
    conn = rmq_conn.open()
    ch = conn.channel()
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
        return 1

    # Queue binding for getting task request from JAWS
    ch.queue_bind(exchange=JGI_JTM_MAIN_EXCH,
                  queue=JTM_TASK_REQUEST_Q,
                  routing_key=JTM_TASK_REQUEST_Q)

    logger.info("JGI Task Manager, version: %s" % (VERSION))
    logger.info("JTM config file: %s" % (CONFIG.config_file))
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

    # MySQL: prepare task table
    db = DbSqlMysql(db=MYSQL_DB, config=CONFIG)
    db.ddl(JTM_SQL["create_database"] % MYSQL_DB)
    db.ddl(JTM_SQL["use_database"] % MYSQL_DB)
    db.ddl(JTM_SQL["create_table_tasks"])
    db.ddl(JTM_SQL["create_table_runs"])
    db.ddl(JTM_SQL["create_table_workers"])
    db.close()

    logger.debug("Main pid = {}".format(PARENT_PROCESS_ID))

    plist = list()

    # Start heartbeat receiving proc
    worker_hb_queue_name = WORKER_HB_Q_POSTFIX
    try:
        rect_hb_from_worker_proc_hdl = mp.Process(target=recv_hb_from_worker_proc,
                                                  args=(worker_hb_queue_name,
                                                        log_dir_name,
                                                        b_resource_usage_log_on))
        rect_hb_from_worker_proc_hdl.start()
        plist.append(rect_hb_from_worker_proc_hdl)
    except Exception as e:
        logger.exception("recv_hb_from_worker_proc: {}".format(e))
        proc_clean(plist)
        conn_clean(conn, ch)
        sys.exit(1)

    logger.info("Waiting for worker\'s heartbeats from %s", worker_hb_queue_name)

    # Start cancelled task cleanup proc
    try:
        task_kill_proc_hdl = mp.Process(target=task_kill_proc)
        task_kill_proc_hdl.start()
        plist.append(task_kill_proc_hdl)
    except Exception as e:
        logger.exception("recv_result_from_workers_proc: {}".format(e))
        proc_clean(plist)
        conn_clean(conn, ch)
        sys.exit(1)

    # Start worker cleanup proc
    try:
        worker_cleanup_proc_hdl = mp.Process(target=zombie_worker_cleanup_proc)
        worker_cleanup_proc_hdl.start()
        plist.append(worker_cleanup_proc_hdl)
    except Exception as e:
        logger.exception("recv_result_from_workers_proc: {}".format(e))
        proc_clean(plist)
        conn_clean(conn, ch)
        sys.exit(1)

    # Start result receiving proc
    # listening to JTM_INNER_RESULT_Q to which all workers will send result messages
    try:
        recv_result_from_worker_proc_hdl = mp.Process(target=recv_result_from_workers_proc)
        recv_result_from_worker_proc_hdl.start()
        plist.append(recv_result_from_worker_proc_hdl)
    except Exception as e:
        logger.exception("recv_result_from_workers_proc: {}".format(e))
        proc_clean(plist)
        conn_clean(conn, ch)
        sys.exit(1)

    # Checking the total number of child processes
    try:
        check_processes_hdl = mp.Process(target=check_processes,
                                         args=(plist,))
        check_processes_hdl.start()
        plist.append(check_processes_hdl)
    except Exception as e:
        logger.exception("check_processes: {}".format(e))
        proc_clean(plist)
        conn_clean(conn, ch)
        sys.exit(1)

    def signal_handler(signum, frame):
        proc_clean(plist)

    signal.signal(signal.SIGTERM, signal_handler)

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
        proc_clean(plist)
        conn_clean(conn, ch)
        sys.exit(1)

    # Keep consuming messages from task request queue
    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        ch.stop_consuming()

    # unreachable
    if conn:
        conn.close()

    return 0
