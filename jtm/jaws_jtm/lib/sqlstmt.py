#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)


JTM_SQL = {
    "create_database": """
        CREATE DATABASE IF NOT EXISTS %s;""",
    "set_timezone": """
        SET GLOBAL time_zone = '-7:00';""",
    "use_database": """
        USE %s;""",
    "create_table_tasks": """
        CREATE TABLE if not exists tasks(
            taskId INTEGER PRIMARY KEY AUTO_INCREMENT,
            userCmd VARCHAR(1024) NOT NULL,
            outFiles VARCHAR(1024) NOT NULL,
            doneFlag VARCHAR(2) NOT NULL DEFAULT '0',
            retryCnt INTEGER NOT NULL DEFAULT 0,
            taskType TINYINT NOT NULL
        ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
    "create_table_runs": """
        CREATE TABLE if not exists runs(
            runId INTEGER PRIMARY KEY AUTO_INCREMENT,
            taskId INTEGER NOT NULL,
            status INTEGER NOT NULL,
            resource VARCHAR(1024),
            workerId2 INTEGER NOT NULL DEFAULT 0,
            childPid INTEGER,
            cancelled TINYINT NOT NULL DEFAULT 0,
            queueDate TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
            startDate TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
            endDate TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (taskId) REFERENCES tasks(taskId) ON DELETE CASCADE,
            FOREIGN KEY (taskId) REFERENCES tasks(taskId) ON UPDATE CASCADE,
            UNIQUE KEY taskid_idx (taskId),
            KEY workerid2_idx (workerId2)
        ) ENGINE= InnoDB DEFAULT CHARSET=latin1;""",
    "create_table_workers": """
        CREATE TABLE if not exists workers(
            workerId2 INTEGER PRIMARY KEY AUTO_INCREMENT,
            workerId CHAR(24),
            workerType TINYINT,
            slurmJobId INTEGER NOT NULL DEFAULT 0,
            startDate TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
            endDate TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
            lifeLeft SMALLINT NOT NULL DEFAULT 0,
            memPerNode VARCHAR(10),
            memPerCore VARCHAR(10),
            nCores TINYINT,
            jobTime VARCHAR(8),
            cloneRate FLOAT,
            cloneCnt TINYINT NOT NULL DEFAULT 0,
            hostName VARCHAR(255),
            jtmHostName VARCHAR(255),
            ipAddress VARCHAR(15),
            poolName VARCHAR(255),
            nWorkersPerNode TINYINT,
            KEY slurmjid_idx (slurmJobId),
            UNIQUE KEY workerkey_idx (workerId)
        ) ENGINE=InnoDB DEFAULT CHARSET=latin1;""",
    "select_exists_workers_by_workerid": """
        SELECT EXISTS (
        SELECT 1
        FROM workers
        WHERE workerId='%(worker_id)s');""",
    "insert_workers_workerid_slurmjobid": """
        INSERT INTO workers (workerId,
                             slurmJobId,
                             workerType,
                             nWorkersPerNode)
        VALUES ('%(worker_id)s',
                %(slurm_jid)d,
                %(worker_type)d,
                %(nwpn)d);""",
    "insert_workers_workerid_workertype_poolname": """
        INSERT INTO workers (workerId,
                             workerType,
                             poolName,
                             jtmHostName,
                             lifeLeft,
                             slurmJobId,
                             nWorkersPerNode)
        VALUES ('%(worker_id)s',
            %(worker_type)d,
            '%(pool_name)s',
            '%(jtm_host_name)s',
            %(lifeleft)d,
            %(slurm_jid)d,
            %(nwpn)d);""",
    "update_workers_enddate_lifeleft_by_workerid": """
        UPDATE workers
        SET endDate ='%(now)s',
            lifeLeft = %(life_left)d,
            memPerNode = '%(mpn)s',
            memPerCore = '%(mpc)s',
            nCores = %(num_cores)d,
            jobTime = '%(job_time)s',
            cloneRate = %(clone_rate)f,
            hostName = '%(host_name)s',
            jtmHostName = '%(jtm_host_name)s',
            ipAddress = '%(ipaddr)s',
            poolName = '%(pool_name)s',
            slurmJobId = %(slurm_jid)d
        WHERE workerId = '%(worker_id)s';""",
    "update_runs_status_by_taskid": """
        UPDATE runs
        SET status = %(status_id)d,
            workerId2 = (select workerId2
                         from workers
                         where workerId='%(worker_id)s'),
            childPid = %(child_pid)d,
            resource = '%(resource_log)s'
        WHERE taskId = %(task_id)d
              and status < 2;""",
    "select_clonecnt_workers_by_workerid": """
        SELECT cloneCnt
        FROM workers
        WHERE workerId = '%(worker_id)s';""",
    "update_workers_clonecnt_by_slurmjobid": """
        UPDATE workers
        set cloneCnt = cloneCnt + 1
        WHERE slurmJobId = '%(slurm_jid)d';""",
    "select_workerid_workers_by_lifeleft_jtmhostname": """
        SELECT workerId
        FROM workers
        WHERE lifeLeft >= 0
              and jtmHostName = '%(jtm_host_name)s'""",
    "select_sum_nwpn_workers_by_lifeleftt_jtmhostname": """
        SELECT sum(nWorkersPerNode)
        FROM workers
        WHERE lifeLeft >= 0
              and jtmHostName = '%(jtm_host_name)s'""",
    "update_workers_lifeleft_by_workerid": """
        UPDATE workers
        SET lifeLeft = -1
        WHERE workerId = '%(worker_id)s'
              and jtmHostName = '%(jtm_host_name)s'""",
    "select_chilpid_runs_by_tid": """
        SELECT childPid
        FROM runs
        WHERE taskId = %(task_id)d""",
    "update_workers_lifeleft_for_last": """
        UPDATE workers
        SET lifeLeft = -1
        WHERE lifeLeft >= 0""",
    "select_workerid2_workers_by_wid": """
        select workerId2
        from workers
        where workerId = '%(worker_id)s'""",
    "update_runs_status_workerid2_by_taskid_2": """
        UPDATE runs
        SET status = %(status_id)d,
            workerId2 = %(wid2)d,
            endDate = '%(now)s'
        WHERE taskId = %(task_id)d;""",
    # "select_status_runs_by_taskid": """
    #     SELECT status
    #     FROM runs
    #     WHERE taskId = %(task_id)d""",
    "select_resource_runs_by_taskid": """
        SELECT resource
        FROM runs
        WHERE taskId = %(task_id)d""",
    "select_workerid_workers_by_tid": """
        SELECT workerId
        FROM workers as w
        INNER JOIN runs as r on r.workerId2 = w.workerId2
        WHERE taskId = %(task_id)d;""",
    "update_tasks_doneflag_by_taskid": """
        UPDATE tasks
        SET doneFlag = %(done_flag)s
        WHERE taskId = %(task_id)s;""",
    "insert_tasks_usercmd_outfiles": """
        INSERT INTO tasks (userCmd, outFiles, doneFlag, retryCnt, taskType)
        VALUES ('%s', '%s', '%s', %d, %d);""",
    "select_last_insert_id": """
        SELECT LAST_INSERT_ID();""",
    "insert_runs_tid_sid": """
        INSERT INTO runs (taskId, status)
        VALUES (%(task_id)d, %(status_id)d);""",
    "update_runs_tid_startdate_by_tid": """
        UPDATE runs
        SET status = %(status_id)d, startdate = '%(now)s'
        WHERE taskId = %(task_id)d
              and status < 1;""",
    "update_runs_status_to_terminated_by_tid": """
        UPDATE runs
        SET status = %(status_id)d, enddate = '%(now)s'
        WHERE taskId = %(task_id)d;""",
    "update_runs_status_cancelled_to_terminated_by_tid": """
        UPDATE runs
        SET status = %(status_id)d, enddate = '%(now)s', cancelled = 2
        WHERE taskId = %(task_id)d;""",
    "update_runs_cancelled_by_tid": """
        UPDATE runs
        SET cancelled = 1, enddate = '%(now)s'
        WHERE taskId = %(task_id)d;""",
    "select_status_runs_by_taskid": """
        SELECT status
        FROM runs
        WHERE taskId = %(task_id)d;""",
    "select_status_cancelled_runs_by_taskid": """
        SELECT status, cancelled
        FROM runs
        WHERE taskId = %(task_id)d;""",
    "select_cancelled_runs_by_taskid": """
        SELECT cancelled
        FROM runs
        WHERE taskId = %(task_id)d;""",
    "select_tids_runs_by_cancelled_and_wid": """
        SELECT taskId
        FROM runs
        WHERE cancelled = 1
              AND workerId2 != 0
              AND childPid != 0
              AND status != -4;""",
    "select_tids_runs_by_cancelled_and_wid2": """
        SELECT taskId
        FROM runs
        WHERE cancelled = 1
              AND workerId2 != 0
              AND childPid != 0;""",
    "select_count_workers_by_poolname_enddate": """
        SELECT count(*)
        FROM workers
        WHERE poolName = '%(pool_name)s'
              and ((lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
                  or
                  (slurmJobId > 1 and startDate = endDate));""",
    "select_count_workers_by_jtm_host_name": """
    SELECT count(*)
    FROM workers
    WHERE jtmHostName = '%(jtm_host_name)s'
          and ((lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
              or
              (slurmJobId > 1 and startDate = endDate));""",
    "select_count_workers_by_jtm_host_name_poolname_enddate": """
    SELECT count(*)
    FROM workers
    WHERE jtmHostName = '%(jtm_host_name)s'
          and poolName = '%(pool_name)s'
          and ((lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
              or
              (slurmJobId > 1 and startDate = endDate));""",
    "select_sum_nwpn_workers_by_poolname_enddate": """
    SELECT SUM(nWorkersPerNode)
    FROM workers
    WHERE poolName = '%(pool_name)s'
          and (lifeLeft >= 0
              or
              (lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
              or
              (slurmJobId > 1 and startDate = endDate));""",
    "select_sum_nwpn_workers_by_jtm_host_name_enddate": """
    SELECT SUM(nWorkersPerNode)
    FROM workers
    WHERE jtmHostName = '%(jtm_host_name)s'
          and (lifeLeft >= 0
              or
              (lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
              or
              (slurmJobId > 1 and startDate = endDate));""",
    "select_sum_nwpn_workers_by_jtm_host_name_poolname_enddate": """
        SELECT SUM(nWorkersPerNode)
        FROM workers
        WHERE jtmHostName = '%(jtm_host_name)s'
              and poolName = '%(pool_name)s'
              and (lifeLeft >= 0
                  or
                  (lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
                  or
                  (slurmJobId > 1 and startDate = endDate));""",
    "select_count_distinct_jid_workers_by_poolname": """
        SELECT COUNT(DISTINCT(slurmJobId))
        FROM workers
        WHERE poolName = '%(pool_name)s'
              and workerType = 2
              and ((lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < 18)
                  or
                  ((slurmJobId > 1 and startDate = endDate)
                  or
                  (slurmJobId > 1
                   and lifeLeft = -1
                   and TIMESTAMPDIFF(SECOND, endDate, NOW()) < 18)
                   and TIMESTAMPDIFF(SECOND, endDate, NOW()) >= 0)
        );""",
    "select_distinct_jid_workers_by_poolname": """
        SELECT DISTINCT(slurmJobId)
        FROM workers
        WHERE poolName = '%(pool_name)s'
              and workerType = 2
              and ((lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < 18)
                  or
                  ((slurmJobId > 1 and startDate = endDate)
                  or
                  (slurmJobId > 1
                   and lifeLeft = -1
                   and TIMESTAMPDIFF(SECOND, endDate, NOW()) < 18)
                   and TIMESTAMPDIFF(SECOND, endDate, NOW()) >= 0)
        );""",
    # process_remove_pool()
    "update_lifeleft_enddate_workers_by_poolname": """
        update workers
        set lifeLeft = -1, endDate = '%(now)s'
        where lifeLeft != -1
              and poolName = '%(pool_name)s'
              and workerId2 > 0;""",
    "select_jid_workers_by_poolname": """
        SELECT DISTINCT(slurmJobId)
        FROM workers
        where workerType = 2
              and slurmJobId > 1
              and (lifeLeft > 0 or lifeLeft = -2)
              and poolName = '%(pool_name)s';""",
    "select_all_jid_workers_by_poolname": """
        SELECT DISTINCT(slurmJobId)
        FROM workers
        where slurmJobId > 1
              and poolName = '%(pool_name)s';""",
    "select_slurmjid_workers_by_lifeleft": """
        select DISTINCT(slurmJobId)
        from workers
        where lifeLeft = -2
              and slurmJobId > 0;""",
    # zombie_worker_cleanup_thead()
    "update_workers_lifeleft_by_slurmjid": """
        update workers
        set lifeLeft = -1, endDate = '%(now)s'
        where slurmJobId = %(slurm_jid)d;""",
    "delete_from_workers_by_poolname": """
        delete from workers
        where poolName = '%(pool_name)s';""",
    "delete_from_workers_by_slurmjid": """
        delete from workers
        where slurmJobId = %(slurm_jid)d;""",
    "update_runs_status_cancelled_by_status_workerId2_poolname": """
        update runs
        set status = -4, cancelled = 3
        where status in (0, 1, 2) and
              workerId2 in (select workerId2
                            from workers
                            where poolName = '%(pool_name)s';""",
    "update_runs_status_cancelled_by_status_workerId2_jid": """
        update runs
        set status = -4, cancelled = 4
        where status in (0, 1, 2) and
              workerId2 in (select workerId2
                            from workers
                            where slurmJobId = %(slurm_jid)d;"""
}
