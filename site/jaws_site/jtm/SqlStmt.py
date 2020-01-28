#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Copyright 2018-2019 (C) DOE JGI, LBL
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
        WHERE workerId='%(wid)s');""",

    "insert_workers_workerid_slurmjobid": """
        INSERT INTO workers (workerId, slurmJobId, workerType, nWorkersPerNode)
        VALUES ('%(wid)s', %(slurmjid)d, %(wtype)d, %(nwpn)d);""",

    "insert_workers_workerid_workertype_poolname": """
        INSERT INTO workers (workerId, workerType, poolName, jtmHostName, lifeLeft, slurmJobId, nWorkersPerNode)
        VALUES ('%(wid)s', %(wtype)d, '%(poolname)s', '%(jtmhostname)s', %(lifeleft)d, %(slurmjobid)d, %(nwpn)d);""",

    "update_workers_enddate_lifeleft_by_workerid": """
        UPDATE workers
        SET endDate='%(now)s',
            lifeLeft=%(left)d,
            memPerNode='%(mpn)s',
            memPerCore='%(mpc)s',
            nCores=%(ncores)d,
            jobTime='%(jtime)s',
            cloneRate=%(clonerate)f,
            hostName='%(hname)s',
            jtmHostName='%(jtmhname)s',
            ipAddress='%(ipaddr)s',
            poolName='%(poolname)s', 
            slurmJobId=%(slurmjid)d
        WHERE workerId='%(wid)s';""",

    "update_runs_status_by_taskid": """
        UPDATE runs
        SET status=%(sid)d, workerId2=(select workerId2 from workers where workerId='%(wid)s'), childPid=%(cpid)d, resource='%(resourcelog)s'
        WHERE taskId=%(tid)d and status < 2;""",

    "select_clonecnt_workers_by_workerid": """
        SELECT cloneCnt
        FROM workers
        WHERE workerId='%(wid)s';""",

    # "update_workers_clonecnt_by_workerid": """
    #     UPDATE workers
    #     set cloneCnt = cloneCnt + 1
    #     WHERE workerId='%(wid)s';""",

    "update_workers_clonecnt_by_slurmjobid": """
        UPDATE workers
        set cloneCnt = cloneCnt + 1
        WHERE slurmJobId='%(slurmjobid)d';""",

    # "select_workerid_workers_by_lifeleft": """
    #     SELECT workerId
    #     FROM workers
    #     WHERE lifeLeft >= 0""",

    "select_workerid_workers_by_lifeleft_jtmhostname": """
        SELECT workerId
        FROM workers
        WHERE lifeLeft >= 0 and jtmHostName='%(jtmhostname)s'""",

    # "select_sum_nwpn_workers_by_lifeleft": """
    #     SELECT sum(nWorkersPerNode)
    #     FROM workers
    #     WHERE lifeLeft >= 0""",

    "select_sum_nwpn_workers_by_lifeleftt_jtmhostname": """
        SELECT sum(nWorkersPerNode)
        FROM workers
        WHERE lifeLeft >= 0 and jtmHostName='%(jtmhostname)s'""",

    "update_workers_lifeleft_by_workerid": """
        UPDATE workers
        SET lifeLeft=-1
        WHERE workerId='%(wid)s' and jtmHostName='%(jtmhostname)s'""",

    "select_chilpid_runs_by_tid": """
        SELECT childPid
        FROM runs
        WHERE taskId = %(tid)d""",

    "update_workers_lifeleft_for_last": """
        UPDATE workers
        SET lifeLeft=-1
        WHERE lifeLeft>=0""",

    "select_workerid2_workers_by_wid": """
        select workerId2 
        from workers 
        where workerId='%(wid)s'""",

    # "update_runs_status_workerid2_by_taskid": """
    #     UPDATE runs
    #     SET status=%(sid)d, workerId2=(SELECT workerId2 FROM workers WHERE workerId='%(wid)s'), endDate='%(now)s'
    #     WHERE taskId=%(tid)d;""",

    "update_runs_status_workerid2_by_taskid_2": """
        UPDATE runs
        SET status=%(sid)d, workerId2=%(wid2)d, endDate='%(now)s'
        WHERE taskId=%(tid)d;""",

    "select_status_runs_by_taskid": """
        SELECT status
        FROM runs
        WHERE taskId = %(tid)d""",

    "select_resource_runs_by_taskid": """
        SELECT resource
        FROM runs
        WHERE taskId = %(tid)d""",

    "select_workerid_workers_by_tid": """
        SELECT workerId
        FROM workers as w
        INNER JOIN runs as r on r.workerId2 = w.workerId2
        WHERE taskId = %(tid)d;""",

    "update_tasks_doneflag_by_taskid": """
        UPDATE tasks
        SET doneFlag=%(dflag)s
        WHERE taskId=%(tid)s;""",

    "insert_tasks_usercmd_outfiles": """
        INSERT INTO tasks (userCmd, outFiles, doneFlag, retryCnt, taskType)
        VALUES ('%s', '%s', '%s', %d, %d);""",

    "select_last_insert_id": """
        SELECT LAST_INSERT_ID();""",

    "insert_runs_tid_sid": """
        INSERT INTO runs (taskId, status)
        VALUES (%(tid)d, %(sid)d);""",

    "update_runs_tid_startdate_by_tid": """
        UPDATE runs
        SET status=%(sid)d, startdate='%(now)s'
        WHERE taskId=%(tid)d and status < 1;""",

    "update_runs_status_to_terminated_by_tid": """
        UPDATE runs
        SET status=%(sid)d, enddate='%(now)s'
        WHERE taskId=%(tid)d;""",

    "update_runs_cancelled_by_tid": """
        UPDATE runs
        SET cancelled=1, enddate='%(now)s'
        WHERE taskId=%(tid)d;""",

    "select_status_runs_by_taskid": """
        SELECT status
        FROM runs
        WHERE taskId=%(tid)d;""",

    "select_cancelled_runs_by_taskid": """
        SELECT cancelled
        FROM runs
        WHERE taskId=%(tid)d;""",

    "select_tids_runs_by_cancelled_and_wid": """
        SELECT taskId
        FROM runs
        WHERE cancelled=1 AND workerId2!=0 AND childPid!=0 AND status!=-4;""",

    # "select_workerid_workers_by_taskid": """
    #     SELECT workerId
    #     FROM workers as w
    #     LEFT JOIN runs as r on r.workerId2 = w.workerId2
    #     WHERE taskId = %(tid)d;""",

    # "select_count_workers_by_poolname_enddate": """
    #     SELECT count(*)
    #     FROM workers
    #     WHERE lifeLeft != -1 and
    #         poolName='%(poolname)s' and
    #         TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d);""",

    "select_count_workers_by_poolname_enddate": """
        SELECT count(*)
        FROM workers
        WHERE poolName='%(poolname)s' and (
                (lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d) 
                or 
                (slurmJobId > 1 and startDate = endDate));""",

    "select_sum_nwpn_workers_by_poolname_enddate": """
        SELECT SUM(nWorkersPerNode)
        FROM workers
        WHERE poolName='%(poolname)s' and (
                (lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d) 
                or 
                (slurmJobId > 1 and startDate = endDate));""",

    # "select_count_workers_by_jtm_host_name": """
    #     SELECT count(*)
    #     FROM workers
    #     WHERE lifeLeft != -1 and
    #         jtmHostName='%(jtmhostname)s' and
    #         TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d;""",

    "select_count_workers_by_jtm_host_name": """
    SELECT count(*)
    FROM workers
    WHERE jtmHostName='%(jtmhostname)s' and (
            (lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
            or 
            (slurmJobId > 1 and startDate = endDate));""",

    "select_sum_nwpn_workers_by_jtm_host_name": """
    SELECT SUM(nWorkersPerNode)
    FROM workers
    WHERE jtmHostName='%(jtmhostname)s' and (
            (lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
            or 
            (slurmJobId > 1 and startDate = endDate));""",
                                         
    # "select_count_workers_by_poolname": """
    #     SELECT count(*)
    #     FROM workers
    #     WHERE lifeLeft != -1 and
    #     poolName='%(poolname)s' and
    #     slurmJobId != 0 and
    #     (startDate = endDate or TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d);"""

    # "select_count_workers_by_poolname": """
    #     SELECT count(*)
    #     FROM workers
    #     WHERE poolName='%(poolname)s' and (
    #         (lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
    #         or
    #         (slurmJobId > 1 and startDate = endDate)
    #         or
    #         (slurmJobId > 1 and lifeLeft = -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
    #     );""",

    "select_count_distinct_jid_workers_by_poolname": """
        SELECT COUNT(DISTINCT(slurmJobId))
        FROM workers
        WHERE poolName='%(poolname)s' and workerType = 2 and (
            (lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < 18)
            or
            ((slurmJobId > 1 and startDate = endDate)
            or
            (slurmJobId > 1 and lifeLeft = -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < 18) and TIMESTAMPDIFF(SECOND, endDate, NOW()) >= 0)
        );""",

    # "select_count_distinct_workerid_workers_by_poolname": """
    #     SELECT COUNT(workerId)
    #     FROM workers
    #     WHERE poolName='%(poolname)s' and (
    #         (lifeLeft != -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
    #         or
    #         (slurmJobId > 1 and startDate = endDate)
    #         or
    #         (slurmJobId > 1 and lifeLeft = -1 and TIMESTAMPDIFF(SECOND, endDate, NOW()) < %(hbinterval)d)
    #     );""",

    # process_remove_pool()
    "update_lifeleft_enddate_workers_by_poolname": """
        update workers 
        set lifeLeft = -1, endDate = '%(now)s'
        where lifeLeft != -1 and poolName = '%(poolname)s' and workerId2 > 0;""",

    "select_jid_workers_by_poolname": """
        SELECT DISTINCT(slurmJobId)
        FROM workers
        where workerType = 2 and slurmJobId > 1 and (lifeLeft > 0 or lifeLeft = -2) and poolName = '%(poolname)s';""",

    "select_slurmjid_workers_by_lifeleft": """
        select DISTINCT(slurmJobId) 
        from workers 
        where lifeLeft = -2 and slurmJobId > 0;""",

    # zombie_worker_cleanup_thead()
    "update_workers_lifeleft_by_slurmjid": """
        update workers 
        set lifeLeft = -1, endDate = '%(now)s'
        where slurmJobId = %(slurmjid)d"""

}