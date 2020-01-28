#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Copyright 2018-2019 (C) DOE JGI, LBL
# Seung-Jin Sul (ssul@lbl.gov)

# import base64
import os
import socket
import pika

PIKA_VER = pika.__version__

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------
VERSION = "5.4.6"
#PRODUCTION = True  # prod
PRODUCTION = False  # dev

################################################################################
# Specify your way of env activation
################################################################################
# This should be set as well in jtm.conf for Cromwell
# If conda,
# ENV_ACTIVATION = "module load python; source activate jtm"
# If conda by jaws,
# ENV_ACTIVATION = "source activate jaws"
# If virtualenv
ENV_ACTIVATION = "source ~/venv/bin/activate" if PRODUCTION else "source ~/venv-dev/bin/activate"

# Accnt name
# This is for running one or more jtm systems under different names
# We cannot  simply use getpass.getuser(), b/c some users have different user names across systems
# Also getpass.getuser() is invalid in aws
#
# import getpass
# USER_NAME = getpass.getuser()
USER_NAME = "jtm" if PRODUCTION else "jtm_dev"

# This is to set the location of jtm instance
# import socket
# HOST_NAME = socket.gethostname()
if "NERSC_HOST" in os.environ:
    JTM_HOST_NAME = os.environ["NERSC_HOST"]
elif "JTM_HOST_NAME" in os.environ:  # for custom name like ["aws' | 'olcf' | 'pc"]
    JTM_HOST_NAME = os.environ["JTM_HOST_NAME"]
elif "HOSTNAME" in os.environ:
    JTM_HOST_NAME = os.environ["HOSTNAME"]
else:
    JTM_HOST_NAME = socket.gethostname().replace(".", "_")

assert JTM_HOST_NAME is not None

# Available cluster/cloud
COMPUTE_RESOURCES = ["denovo", "cori", "lawrencium", "rhea", "ssul-dm.dhcp.lbl.gov", "aws", "jaws_lbl_gov"]

# Note: jtm also uses this account name for ssh-sbatching dynamic workers on denovo/cori
CNAME = JTM_HOST_NAME + "." + USER_NAME
if JTM_HOST_NAME == "cori":
    if PRODUCTION:
        # prod
        JOB_LOG = os.path.join(os.environ["CSCRATCH"], "jtm/%s/jobs" % (CNAME))  # save slurm related logs
        JTM_LOG = os.path.join(os.environ["CSCRATCH"], "jtm/%s/logs" % (CNAME))  # save jtm and jtm-worker related logs
    else:
        # dev
        JOB_LOG = os.path.join(os.environ["CSCRATCH"], "dev/jtm/%s/jobs" % (CNAME))
        JTM_LOG = os.path.join(os.environ["CSCRATCH"], "dev/jtm/%s/logs" % (CNAME))
else:
    if PRODUCTION:
        # prod
        JOB_LOG = os.path.join(os.environ["SCRATCH"], "jtm/%s/jobs" % (CNAME))  # save slurm related logs
        JTM_LOG = os.path.join(os.environ["SCRATCH"], "jtm/%s/logs" % (CNAME))  # save jtm and jtm-worker related logs
    else:
        # dev
        JOB_LOG = os.path.join(os.environ["SCRATCH"], "dev/jtm/%s/jobs" % (CNAME))  # save slurm related logs
        JTM_LOG = os.path.join(os.environ["SCRATCH"], "dev/jtm/%s/logs" % (CNAME))  # save jtm and jtm-worker related logs


################################################################################
# RMQ setup
################################################################################
# RMQ_HOST = "mq.nersc.gov"
# RMQ_USER = "sulsj"
# RMQ_PASS = "synchrotron"
# RMQ_VHOST = "jgi"
# RMQ_PORT = "5672"

# SPIN http://rmq.nersc.gov:60039/
RMQ_HOST = "rmq.nersc.gov"
RMQ_USER = "sulsj"
RMQ_PASS = "telepath-downsize-maxim"
RMQ_VHOST = "jgi"
if JTM_HOST_NAME in ("cori", "denovo"):
    RMQ_PORT = "5672"
else:
    RMQ_PORT = "60042"

# Between jtm and workers
JTM_INNER_MAIN_EXCH = "jgi_jtm_inner_main_exchange"
JTM_WORKER_HB_EXCH = "jgi_jtm_worker_hb_exchange"
JTM_CLIENT_HB_EXCH = "jgi_jtm_client_hb_exchange"
if not PRODUCTION:
    JTM_INNER_MAIN_EXCH = "jgi_jtm_inner_main_exchange_dev"
    JTM_WORKER_HB_EXCH = "jgi_jtm_worker_hb_exchange_dev"
    JTM_CLIENT_HB_EXCH = "jgi_jtm_client_hb_exchange_dev"

JTM_INNER_REQUEST_Q = "_jtm_inner_request_queue" + "." + CNAME
JTM_INNER_RESULT_Q = "_jtm_inner_result_queue" + "." + CNAME
WORKER_HB_Q_POSTFIX = "_jtm_worker_hb_queue" + "." + CNAME  # worker hb queue
CLIENT_HB_Q_POSTFIX = "_jtm_client_hb_queue" + "." + CNAME  # client hb queue

# Between jaws and jtm
JGI_JTM_MAIN_EXCH = "jgi_jtm_main"   # exchange for user task request and replay b/w JWS and JTM
if not PRODUCTION:
    JGI_JTM_MAIN_EXCH = "jgi_jtm_main_dev"

JTM_TASK_REQUEST_Q = "_jtm_task_request_queue" + "." + CNAME
JTM_TASK_RESULT_Q = "_jtm_task_result_queue" + "." + CNAME

# Poison exchange
JTM_WORKER_POISON_EXCH = "jgi_jtm_poison_exchange"
if not PRODUCTION:
    JTM_WORKER_POISON_EXCH = "jgi_jtm_poison_exchange_dev"
JTM_WORKER_POISON_Q = "_jtm_worker_poison" + "." + CNAME

# Task kill
JTM_TASK_KILL_EXCH = "jgi_jtm_task_kill_exchange"
if not PRODUCTION:
    JTM_TASK_KILL_EXCH = "jgi_jtm_task_kill_exchange_dev"
JTM_TASK_KILL_Q = "_jtm_task_kill" + "." + CNAME

################################################################################
# JTM system variables
################################################################################

# In Manager
CLIENT_HB_RECV_INTERVAL = 6.0   # worker->client heartbeat receiving interval
CLIENT_HB_SEND_INTERVAL = 5.0   # client->worker heartbeat sending interval
RESULT_RECEIVE_INTERVAL = 0.1   # interval to check the result from jtm
TASK_KILL_INTERVAL = 5.0        # manager's interval to check if there is task to terminate
WORKER_KILL_INTERVAL = 300.0
WORKER_INFO_UPDATE_WAIT = 0.5   # to ensure live worker info is updated in workers' table
RUNS_INFO_UPDATE_WAIT = 0.5     # to ensure updating runs table
WORKER_HB_CHECK_MAX_COUNT = 0   # live worker checking count limit
TASK_STAT_UPDATE_INTERVAL = 0.5
NUM_RESULT_RECV_THREADS = 5     # num threads to recv results from workers
NUM_TASK_REQUEST_THREADS = 10   # num threads to recv tasks from jtm-*

# In Worker
WORKER_HB_RECV_INTERVAL = 6.0   # client->worker heartbeat receiving interval in worker
WORKER_HB_SEND_INTERVAL = 3.0   # worker->client heartbeat sending interval in worker
WORKER_TIMEOUT = 0              # timeout for waiting the client's heartbeat; 0=no timeout
FILE_CHECKING_MAX_TRIAL = 3     # max number of trial for checking
FILE_CHECK_INTERVAL = 3.0       # sleep time between output file checking before retrial
FILE_CHECK_INT_INC = 1.5        # increase amount of wait time for file checking

# Options
JTM_STREAM_LOGGING = True         # enable/disable stream logging
JTM_WORKER_STREAM_LOGGING = True  # enable/disable stream logging
JTM_FILE_LOGGING = True           # enable/disable file logging
JTM_WORKER_FILE_LOGGING = True    # enable/disable file logging

# JTM Interface
JTMINTERFACE_TIMEOUT = 1200    # jtm_* CLI command timeout in secs.
JTMINTERFACE_MAX_TRIAL = 1200  # jtm_* CLI command timeout in secs.

# Not used
# CLIENT_HB_RECEIVE_INITIAL_INTERVAL = 1  # heartbeat receiving interval in sec of client
# CLIENT_HB_RECEIVE_INT_INC_MUL = 2.0     # heartbeat increase multiplier of client
# CLIENT_HB_RECEIVE_INT_INC_RATE = 1.2    # heartbeat increase rate of client
# MAX_RETRY_CNT_CHILDPID = 5              # max wait time for valid child pid for a user task
# QTTL = 60000                            #  10 secs ttl for queue


################################################################################
# DB setup
################################################################################
# SQLITE_DB = "jtm.db"

#MYSQL_HOST = "localhost"
#MYSQL_PORT = 3306
#MYSQL_USER = "jtm"
#MYSQL_PW = "jtm"

# MYSQL_HOST = "gpdb23.nersc.gov"
# MYSQL_PORT = 3306
# MYSQL_USER = "sulsj"
# #MYSQL_PW = "epleRicket=42"
# MYSQL_PW = "tialWixOn3##"
# MYSQL_DB = "jtm"

# NOTE: Don't forget to set the time zone of MySQL
#  mysql> SET GLOBAL time_zone = '-7:00';
if PRODUCTION:
    MYSQL_HOST = "db.mysql.prod-cattle.stable.spin.nersc.org"
    MYSQL_PORT = 60006
    MYSQL_PW = "bjwHVl275hSl2sjbbzhwiipljwV"
else:
    MYSQL_HOST = "db.mysql.dev-cattle.stable.spin.nersc.org"
    MYSQL_PORT = 60005
    MYSQL_PW = "klsjBl2kjVl9s#k5V"
    
MYSQL_USER = "ssul"
# MYSQL_DB = "jtm" if PRODUCTION else "jtm_dev"
MYSQL_DB = "jtm"
MYSQL_POOLNAME = "jtm"
MYSQL_POOLSIZE = 10  # db connection pool size


################################################################################
# types
################################################################################
TASK_TYPE = {"task": 0,
             "status": 1,
             "kill": 2,
             "resource": 3,
             "term": 9,  # worker termination
             "hb": 88,  # jtm hb to worker
             "check_manager": 98,  # jgi-task-manager check
             "check_worker": 99,  # jtm-worker check
             "remove_pool": 100
             }

TASK_STATUS = {"ready": 0,
               "queued": 1,
               "running": 2,
               "success": 4,
               "outputerror": -1,    # output file(s) not found.
               "failed": -2,
               "outofresource": -3,  # out of mem
               "terminated": -4,     # terminated by user
               "invalidtask": -5     # task definition in the mesasge from jtm_submit is not valid
               }

DONE_FLAGS = {"success": 1,
              "success with correct output file(s)": 2,
              "failed to check output file(s)": -1,
              "failed to run user command": -2,
              "failed with out-of-mem": -3,  # not used
              "failed with user termination": -4
              }

WORKER_TYPE = {"manual": 0,
               "static": 1,
               "dynamic": 2
               }

HB_MSG = {"child_pid": 1,
          "clone_time_rate": 2,
          "cpu_load": 3,
          "end_date": 4,
          "host_name": 5,
          "ip_address": 6,
          "job_time": 7,
          "life_left": 8,
          "mem_per_core": 9,
          "mem_per_node": 10,
          "num_cores": 11,
          "num_tasks": 12,
          "num_workers_on_node": 13,
          "perc_mem_used": 14,
          "pool_name": 15,
          "ret_msg": 16,
          "rmem_usage": 17,
          "root_pid": 18,
          "run_time": 19,
          "slurm_jobid": 20,
          "task_id": 21,
          "vmem_usage": 22,
          "worker_id": 23,
          "worker_type": 24,
          "jtm_host_name": 25,  # hostname on which jtm manager is running ["denovo" | "cori"]
          "nwpn": 26  # num workers per node
          }

################################################################################
# Job scheduling system setting
################################################################################
#==================
# NERSC Cori
#==================
CLUSTER = "cori"        # default cluster to use
NNODES = 1              # number of nodes for the pool
NWORKERS = 1            # number of workers per node
NCPUS = 1               # num cores for small task
JOBTIME = "00:30:00"    # default wall clock time
MEMPERCPU = "1GB"       # memory pool per core (should be use with -c, ie. core based request)
                        # slurm is configured with this `MaxMemPerCPU=1952MB`

MEMPERNODE = "5GB"
CORI_QOS = "genepool_special"   # ["genepool_special", "genepool_shared", "jgi_shared", "jgi_exvivo"]
CORI_ACCNT = "fungalp"          # cori needs -q, -A
# CORI_ACCNT = "m342"             # for general cori nodes and for KNL
CORI_CONSTRAINT = "haswell"     # [haswell, knl, skylake]
CTR = 0.2                       # cloning time rate = life_left / runtime_requested
                                # ex) if -t=5min, fire sbatch at 4min

#==================
# Lab IT Lawrencium
#==================
LAWRENCIUM_PARTITION = "lr2"
#LAWRENCIUM_QOS = "lr_normal"
LAWRENCIUM_QOS = "lr_debug"
LAWRENCIUM_ACCNT = "pc_jaws"

#==================
# OLCF Rhea
#==================

# EOF
