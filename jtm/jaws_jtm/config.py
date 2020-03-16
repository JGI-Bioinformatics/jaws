#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)
import os
import socket
import pika
import getpass
from cryptography.fernet import Fernet

CRYPT_KEY = b"pRmgMa8T0INjEAfksaq2aafzoZXEuwKI7wDe4c1F8AY="
CIPHER_SUITE = Fernet(CRYPT_KEY)

PIKA_VER = pika.__version__

# -------------------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------------------
VERSION = "5.6.8"
# PRODUCTION = True  # prod
PRODUCTION = False  # dev

# Accnt name
# This is for running one or more jtm systems under different names
# We cannot  simply use getpass.getuser(), b/c some users have different user names across systems
# Also getpass.getuser() is invalid in aws
#

LOGIN_NAME = getpass.getuser()

# Note: the .jtm_custom_user_name file can be created for special debugging
#  with different jtm user name. Jtm user name is either jtm or jtm_dev by default
#  Any string in the .jtm_custom_user_name file will be used for Jtm user name
#  which is useful to debug JTM b/c it won't interfere with the production JTM system.
#
customJtmUserNameFile = os.path.expanduser("~/.jtm")
if os.path.isfile(os.path.expanduser(customJtmUserNameFile)):
    # import getpass
    # USER_NAME = getpass.getuser()
    with open(customJtmUserNameFile, "r") as cjunf:
        USER_NAME = cjunf.readline().strip().replace(".", "_")
else:
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
JTM_HOST_NAME = JTM_HOST_NAME.replace(".", "_")

################################################################################
# Specify your way of env activation
################################################################################

# This should be set as well in Cromwell configuration file
if JTM_HOST_NAME == "cori":
    if PRODUCTION:
        if LOGIN_NAME == "jaws" and JTM_HOST_NAME == "cori":
            # For jaws on cori20
            ENV_ACTIVATION = "source /global/cfs/projectdirs/jaws/jtm/venv/bin/activate"
        else:
            # ENV_ACTIVATION = "source /global/project/projectdirs/jaws_jtm/anaconda3/bin/activate /global/project/projectdirs/jaws_jtm/prod"
            ENV_ACTIVATION = "source /global/cfs/projectdirs/jaws_jtm/anaconda3/bin/activate /global/cfs/projectdirs/jaws_jtm/dev/jtm"
    else:
        if LOGIN_NAME == "jaws" and JTM_HOST_NAME == "cori":
            # For jaws on cori20
            ENV_ACTIVATION = (
                "source /global/cfs/projectdirs/jaws/jtm/venv-dev/bin/activate"
            )
        else:
            # ENV_ACTIVATION = "source /global/project/projectdirs/jaws_jtm/anaconda3/bin/activate /global/project/projectdirs/jaws_jtm/dev"
            ENV_ACTIVATION = "source /global/cfs/projectdirs/jaws_jtm/anaconda3/bin/activate /global/cfs/projectdirs/jaws_jtm/prod/jtm"
elif JTM_HOST_NAME == "jaws_lbl_gov":
    if PRODUCTION:
        ENV_ACTIVATION = "source /global/home/groups-sw/lr_jgicloud/anaconda3/bin/activate /global/home/groups-sw/lr_jgicloud/prod/jtm"
    else:
        ENV_ACTIVATION = "source /global/home/groups-sw/lr_jgicloud/anaconda3/bin/activate /global/home/groups-sw/lr_jgicloud/dev/jtm"
elif JTM_HOST_NAME in ("summit", "rhea", "slate"):
    # OLCF configuration
    pass
else:
    ENV_ACTIVATION = (
        "source ~/venv/bin/activate" if PRODUCTION else "source ~/venv-dev/bin/activate"
    )

DEFAULT_POOL = ["small", "medium", "large", "xlarge"]

# Supported cluster
COMPUTE_RESOURCES = [
    "cori",
    "lawrencium",
    "jgi_cloud",
    "jaws_lbl_gov",
    "jgi_cluster",  # lbl it
    "rhea",
    "summit",  # olcf
    "ssul-dm_dhcp_lbl_gov",
    "sjsul-lm-2_local",
    "sjsul-lm.dhcp.lbnl.us",  # testing
    "summit",
    "rhea",
    "slate",  # OLCF
    "aws",  # future support
]

# $ sacctmgr show assoc user=jaws_jtm -p
# Cluster|Account|User|Partition|Share|Priority|GrpJobs|GrpTRES|GrpSubmit|GrpWall|GrpTRESMins|MaxJobs|MaxTRES|MaxTRESPerNode|MaxSubmit|MaxWall|MaxTRESMins|QOS|Def QOS|GrpTRESRunMins|
# escori|m342|jaws_jtm||1||||||||bb/datawarp=52828800M|||||debug_hsw,debug_knl,flex,interactive,jupyter,long,low_knl,overrun,premium,regular_0,regular_1,regular_bigmem,resv,resv_shared,shared,xfer|||
# cori|m342|jaws_jtm||1||||||||bb/datawarp=52828800M|||||debug_hsw,debug_knl,flex,interactive,jupyter,long,low_knl,overrun,premium,regular_0,regular_1,regular_bigmem,resv,resv_shared,shared,xfer|||
QOS_LIST = [
    "genepool_special",
    "genepool_shared",
    "genepool",
    "regular",
    "jgi_shared",
    "jgi_exvivo",
    "condo_jgicloud",
]

# Note: jtm also uses this account name for ssh-sbatching dynamic workers on cori
CNAME = JTM_HOST_NAME + "." + USER_NAME
if JTM_HOST_NAME == "cori":
    if PRODUCTION:
        # prod
        JOB_LOG = os.path.join(
            os.environ["CSCRATCH"], "jtm/%s/jobs" % (CNAME)
        )  # save slurm related logs
        JTM_LOG = os.path.join(
            os.environ["CSCRATCH"], "jtm/%s/logs" % (CNAME)
        )  # save jtm and jtm-worker related logs
    else:
        # dev
        JOB_LOG = os.path.join(os.environ["CSCRATCH"], "dev/jtm/%s/jobs" % (CNAME))
        JTM_LOG = os.path.join(os.environ["CSCRATCH"], "dev/jtm/%s/logs" % (CNAME))
else:
    if PRODUCTION:
        # prod
        JOB_LOG = os.path.join(
            os.environ["SCRATCH"], "jtm/%s/jobs" % (CNAME)
        )  # save slurm related logs
        JTM_LOG = os.path.join(
            os.environ["SCRATCH"], "jtm/%s/logs" % (CNAME)
        )  # save jtm and jtm-worker related logs
    else:
        # dev
        JOB_LOG = os.path.join(
            os.environ["SCRATCH"], "dev/jtm/%s/jobs" % (CNAME)
        )  # save slurm related logs
        JTM_LOG = os.path.join(
            os.environ["SCRATCH"], "dev/jtm/%s/logs" % (CNAME)
        )  # save jtm and jtm-worker related logs


################################################################################
# RMQ setup
################################################################################
if JTM_HOST_NAME in ("cori"):
    # SPIN server http://rmq.nersc.gov:60039/
    RMQ_HOST = "rmq.nersc.gov"
    RMQ_USER = CIPHER_SUITE.decrypt(
        "gAAAAABdrhu9Cu4sWDP7625IsLLLZgKy3mF_5X-_18QQGzSn4H8D5cix5C5oi4m2EXzvJtL0t6JFemooi1xmB_gLIv88JXshtg==".encode()
    )
    RMQ_PASS = CIPHER_SUITE.decrypt(
        "gAAAAABdrhgBIKwNC0DRX6zZMitFesyOjWgeFuuhHVeTHtWiNfKUhuo6YBIje8l2abqUKCpmtXCbCSKIS0u9J3AgHRu5TPAPvG3rdw4s1V33OsIeEdEfshU=".encode()
    )
    RMQ_VHOST = "jgi"
    if JTM_HOST_NAME in ("cori"):
        RMQ_PORT = "5672"
    else:
        RMQ_PORT = "60042"
elif JTM_HOST_NAME in ("lawrencium", "jgi_cloud", "jaws_lbl_gov"):
    # rmq.lbl.gov
    RMQ_HOST = "rmq.lbl.gov"
    RMQ_USER = CIPHER_SUITE.decrypt(
        "gAAAAABd1ssdmSoy6MtbPocGnbJoD82kXYXI2WkX8hAt976bOe4YBQsz8yQJCiwdQG84TgdmMh2J6i0y5-3j2-KHpcuNelohsQ==".encode()
    )
    RMQ_PASS = CIPHER_SUITE.decrypt(
        "gAAAAABd1ssdmSoy6MtbPocGnbJoD82kXYXI2WkX8hAt976bOe4YBQsz8yQJCiwdQG84TgdmMh2J6i0y5-3j2-KHpcuNelohsQ==".encode()
    )
    RMQ_VHOST = "jaws"
    RMQ_PORT = "5672"
elif JTM_HOST_NAME in ("summit", "rhea", "slate"):
    # OLCF configuration
    RMQ_HOST = "bif113.marble.ccs.ornl.gov"
    RMQ_USER = "admin"
    RMQ_PASS = "********"
    RMQ_VHOST = "/"
    RMQ_PORT = "30449"
else:
    # a docker-compose local rmq broker
    # RMQ_HOST = "localhost"
    # RMQ_USER = "jaws"
    # RMQ_PASS = "jaws"
    # RMQ_VHOST = "/"
    # RMQ_PORT = "5672"
    # rmq.lbl.gov
    RMQ_HOST = "rmq.lbl.gov"
    RMQ_USER = CIPHER_SUITE.decrypt(
        "gAAAAABd1ssdmSoy6MtbPocGnbJoD82kXYXI2WkX8hAt976bOe4YBQsz8yQJCiwdQG84TgdmMh2J6i0y5-3j2-KHpcuNelohsQ==".encode()
    )
    RMQ_PASS = CIPHER_SUITE.decrypt(
        "gAAAAABd1ssdmSoy6MtbPocGnbJoD82kXYXI2WkX8hAt976bOe4YBQsz8yQJCiwdQG84TgdmMh2J6i0y5-3j2-KHpcuNelohsQ==".encode()
    )
    RMQ_VHOST = "jaws"
    RMQ_PORT = "5672"

# Between jtm and workers
# JTM_INNER_MAIN_EXCH = "jgi_jtm_inner_main_exchange"
# JTM_WORKER_HB_EXCH = "jgi_jtm_worker_hb_exchange"
# JTM_CLIENT_HB_EXCH = "jgi_jtm_client_hb_exchange"
# if not PRODUCTION:
#     JTM_INNER_MAIN_EXCH = "jgi_jtm_inner_main_exchange_dev"
#     JTM_WORKER_HB_EXCH = "jgi_jtm_worker_hb_exchange_dev"
#     JTM_CLIENT_HB_EXCH = "jgi_jtm_client_hb_exchange_dev"
JTM_INNER_MAIN_EXCH = "jgi_jtm_inner_main_exch_" + CNAME
JTM_WORKER_HB_EXCH = "jgi_jtm_worker_hb_exch_" + CNAME
JTM_CLIENT_HB_EXCH = "jgi_jtm_client_hb_exch_" + CNAME
if not PRODUCTION:
    JTM_INNER_MAIN_EXCH = "jgi_jtm_inner_main_exch_dev_" + CNAME
    JTM_WORKER_HB_EXCH = "jgi_jtm_worker_hb_exch_dev_" + CNAME
    JTM_CLIENT_HB_EXCH = "jgi_jtm_client_hb_exch_dev_" + CNAME

JTM_INNER_REQUEST_Q = "_jtm_inner_request_queue" + "." + CNAME
JTM_INNER_RESULT_Q = "_jtm_inner_result_queue" + "." + CNAME
WORKER_HB_Q_POSTFIX = "_jtm_worker_hb_queue" + "." + CNAME  # worker hb queue
CLIENT_HB_Q_POSTFIX = "_jtm_client_hb_queue" + "." + CNAME  # client hb queue

# Between jaws and jtm
# JGI_JTM_MAIN_EXCH = "jgi_jtm_main"   # exchange for user task request and replay b/w JWS and JTM
# if not PRODUCTION:
#     JGI_JTM_MAIN_EXCH = "jgi_jtm_main_dev"
JGI_JTM_MAIN_EXCH = (
    "jgi_jtm_main_exch_" + CNAME
)  # exchange for user task request and replay b/w JWS and JTM
if not PRODUCTION:
    JGI_JTM_MAIN_EXCH = "jgi_jtm_main_exch_dev_" + CNAME

JTM_TASK_REQUEST_Q = "_jtm_task_request_queue" + "." + CNAME
JTM_TASK_RESULT_Q = "_jtm_task_result_queue" + "." + CNAME

# Poison exchange
# JTM_WORKER_POISON_EXCH = "jgi_jtm_poison_exchange"
# if not PRODUCTION:
#     JTM_WORKER_POISON_EXCH = "jgi_jtm_poison_exchange_dev"
JTM_WORKER_POISON_EXCH = "jgi_jtm_poison_exch_" + CNAME
if not PRODUCTION:
    JTM_WORKER_POISON_EXCH = "jgi_jtm_poison_exch_dev_" + CNAME
JTM_WORKER_POISON_Q = "_jtm_worker_poison" + "." + CNAME

# Task kill
# JTM_TASK_KILL_EXCH = "jgi_jtm_task_kill_exchange"
# if not PRODUCTION:
#     JTM_TASK_KILL_EXCH = "jgi_jtm_task_kill_exchange_dev"
JTM_TASK_KILL_EXCH = "jgi_jtm_task_kill_exch_" + CNAME
if not PRODUCTION:
    JTM_TASK_KILL_EXCH = "jgi_jtm_task_kill_exch_dev_" + CNAME
JTM_TASK_KILL_Q = "_jtm_task_kill" + "." + CNAME

################################################################################
# JTM system variables
################################################################################

# In Manager
CLIENT_HB_RECV_INTERVAL = 6.0  # worker->client heartbeat receiving interval
# CLIENT_HB_SEND_INTERVAL = 5.0   # client->worker heartbeat sending interval
CLIENT_HB_SEND_INTERVAL = 3.0  # client->worker heartbeat sending interval
RESULT_RECEIVE_INTERVAL = 0.1  # interval to check the result from jtm
TASK_KILL_INTERVAL = 5.0  # manager's interval to check if there is task to terminate
WORKER_KILL_INTERVAL = 300.0  # zombie worker checking interval
WORKER_INFO_UPDATE_WAIT = 0.5  # to ensure live worker info is updated in workers' table
RUNS_INFO_UPDATE_WAIT = 0.5  # to ensure updating runs table
WORKER_HB_CHECK_MAX_COUNT = 0  # live worker checking count limit
TASK_STAT_UPDATE_INTERVAL = 0.5
NUM_RESULT_RECV_THREADS = 5  # num threads to recv results from workers
NUM_TASK_REQUEST_THREADS = 10  # num threads to recv tasks from jtm-*

# In Worker
WORKER_HB_RECV_INTERVAL = (
    6.0  # client->worker heartbeat receiving interval in worker NOT USED!!
)
# WORKER_HB_SEND_INTERVAL = 3.0   # worker->client heartbeat sending interval in worker
WORKER_HB_SEND_INTERVAL = 2.0  # worker->client heartbeat sending interval in worker
WORKER_TIMEOUT = 0  # timeout for waiting the client's heartbeat; 0=no timeout
FILE_CHECKING_MAX_TRIAL = 3  # max number of trial for checking
FILE_CHECK_INTERVAL = 3.0  # sleep time between output file checking before retrial
FILE_CHECK_INT_INC = 1.5  # increase amount of wait time for file checking

# Options
JTM_STREAM_LOGGING = True  # enable/disable stream logging
JTM_WORKER_STREAM_LOGGING = True  # enable/disable stream logging
JTM_FILE_LOGGING = True  # enable/disable file logging
JTM_WORKER_FILE_LOGGING = True  # enable/disable file logging

# JTM Interface
# JTMINTERFACE_TIMEOUT = 300    # jtm_* CLI command timeout in secs. Not used.
JTMINTERFACE_MAX_TRIAL = 300  # jtm_* CLI command timeout in secs.

# Not used
# CLIENT_HB_RECEIVE_INITIAL_INTERVAL = 1  # heartbeat receiving interval in sec of client
# CLIENT_HB_RECEIVE_INT_INC_MUL = 2.0     # heartbeat increase multiplier of client
# CLIENT_HB_RECEIVE_INT_INC_RATE = 1.2    # heartbeat increase rate of client
# MAX_RETRY_CNT_CHILDPID = 5              # max wait time for valid child pid for a user task
# QTTL = 60000                            #  10 secs ttl for queue

JTM_MANAGER_NUM_THREADS = 6  # not used
JTM_WORKER_NUM_THREADS = (
    5  # used in resourceusage.py for getting the number of workers by ps
)

################################################################################
# DB server setup
################################################################################
if JTM_HOST_NAME in ("cori"):
    MYSQL_USER = "ssul"
    MYSQL_DB = "jtm"
    if PRODUCTION:
        MYSQL_HOST = "db.mysql.prod-cattle.stable.spin.nersc.org"
        MYSQL_PORT = 60006
        MYSQL_PW = CIPHER_SUITE.decrypt(
            "gAAAAABdrhip0dSjHxtw6KVgIKTOFWrdAS-TzTqv60WvYN7bMxls5cHaFviIEwfeQByeAsHKQVv6U1oz4gccuwRPfV983rsSAJsuU5YKFIxqJ5XbtDu30so=".encode()
        )
    else:
        MYSQL_HOST = "db.mysql.dev-cattle.stable.spin.nersc.org"
        MYSQL_PORT = 60005
        MYSQL_PW = CIPHER_SUITE.decrypt(
            "gAAAAABdrhoByzbVJVCTH3gV7Wb8aqowrywBx3mccnKvb7nWO-HsMq1nBdteaXW816UIjZBzkWze7ss1vuVEqjVVEAHqsfRJBR_Jp9L1wE9PxZurVNH-1F0=".encode()
        )
elif JTM_HOST_NAME in ("lawrencium", "jgi_cloud", "jaws_lbl_gov"):
    if PRODUCTION:
        MYSQL_USER = "jtm"
        MYSQL_DB = "jtm"
        MYSQL_HOST = "jaws-db.lbl.gov"
        MYSQL_PORT = 3306
        # MYSQL_PW = CIPHER_SUITE.decrypt("gAAAAABd8BAahvN_AYIFL61I0fPCcO0h1Ytilj0EdaO6JPnjdyZBxJqJAIoYTNjHKGDcYjqJB-AnXZcu0gW1Djp0vnemPNTwlirc6YmnLVNnDE380NLzUA4=")
        MYSQL_PW = "*uEF!n!Nb0XT784e^Tf"
    else:
        MYSQL_USER = "jtm_dev"
        MYSQL_DB = "jtm_dev"
        MYSQL_HOST = "jaws-db.lbl.gov"
        MYSQL_PORT = 3306
        # MYSQL_PW = CIPHER_SUITE.decrypt("gAAAAABd8BBNY-vVlNoy_MUt01L3B2gTWoB6n44BlkeB8qQKyC91I1kZr7Tj4j-OSiXzYBG1I1M-knWBN_jS0hZ3JrqTmy1p9NwHLMKauwJV4ug4FWzdbiw=")
        MYSQL_PW = "&g1kT21QhWAp&ZByC9j"
elif JTM_HOST_NAME in ("summit", "rhea", "slate"):
    # Todo: Need OLCF configuration
    pass
else:
    # "ssul-dm_dhcp_lbl_gov", "sjsul-lm-2_local"
    # or other test machines
    if PRODUCTION:
        MYSQL_USER = "jtm"
        MYSQL_DB = "jtm"
        MYSQL_HOST = "jaws-db.lbl.gov"
        MYSQL_PORT = 3306
        # MYSQL_PW = CIPHER_SUITE.decrypt("gAAAAABd8BAahvN_AYIFL61I0fPCcO0h1Ytilj0EdaO6JPnjdyZBxJqJAIoYTNjHKGDcYjqJB-AnXZcu0gW1Djp0vnemPNTwlirc6YmnLVNnDE380NLzUA4=")
        MYSQL_PW = "*uEF!n!Nb0XT784e^Tf"
    else:
        # MYSQL_USER = "jtm_dev"
        # MYSQL_DB = "jtm_dev"
        # MYSQL_HOST = "jaws-db.lbl.gov"
        # MYSQL_PORT = 3306
        # # MYSQL_PW = CIPHER_SUITE.decrypt("gAAAAABd8BBNY-vVlNoy_MUt01L3B2gTWoB6n44BlkeB8qQKyC91I1kZr7Tj4j-OSiXzYBG1I1M-knWBN_jS0hZ3JrqTmy1p9NwHLMKauwJV4ug4FWzdbiw=")
        # MYSQL_PW = "&g1kT21QhWAp&ZByC9j"
        # my local db
        MYSQL_USER = "jtm"
        MYSQL_DB = "jtm"
        MYSQL_HOST = "localhost"
        MYSQL_PORT = 3306
        # MYSQL_PW = CIPHER_SUITE.decrypt("gAAAAABd8BBNY-vVlNoy_MUt01L3B2gTWoB6n44BlkeB8qQKyC91I1kZr7Tj4j-OSiXzYBG1I1M-knWBN_jS0hZ3JrqTmy1p9NwHLMKauwJV4ug4FWzdbiw=")
        MYSQL_PW = "jtm"

# MYSQL_USER = "ssul"
# MYSQL_DB = "jtm"
# MYSQL_POOLNAME = "jtm"
# MYSQL_POOLSIZE = 10  # db connection pool size


################################################################################
# types
################################################################################
TASK_TYPE = {
    "task": 0,
    "status": 1,
    "kill": 2,
    "resource": 3,
    "term": 9,  # worker termination
    "hb": 88,  # jtm hb to worker
    "check_manager": 98,  # jgi-task-manager check
    "check_worker": 99,  # jtm-worker check
    "remove_pool": 100,
}

TASK_STATUS = {
    "ready": 0,
    "queued": 1,
    "running": 2,
    "success": 4,
    "outputerror": -1,  # output file(s) not found.
    "failed": -2,
    "outofresource": -3,  # out of mem
    "terminated": -4,  # terminated by user
    "invalidtask": -5,  # task definition in the mesasge from jtm_submit is not valid
}

DONE_FLAGS = {
    "success": 1,
    "success with correct output file(s)": 2,
    "failed to check output file(s)": -1,
    "failed to run user command": -2,
    "failed with out-of-mem": -3,  # not used
    "failed with user termination": -4,
}

WORKER_TYPE = {"manual": 0, "static": 1, "dynamic": 2}

HB_MSG = {
    "child_pid": 1,
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
    "jtm_host_name": 25,  # hostname on which jtm manager is running ["cori"]
    "nwpn": 26,  # num workers per node, NOT USED
}

################################################################################
# Job scheduling system setting
################################################################################
# ==================
# NERSC Cori
# ==================
CLUSTER = "cori"  # default cluster to use
NNODES = 1  # number of nodes for the pool
NWORKERS = 1  # number of workers per node
NCPUS = 1  # num cores for small task
JOBTIME = "00:30:00"  # default wall clock time
MEMPERCPU = (
    "1GB"  # memory pool per core (should be use with -c, ie. core based request)
)
# slurm is configured with this `MaxMemPerCPU=1952MB`

MEMPERNODE = "5GB"
CORI_QOS = "genepool_special"  # ["genepool_special", "genepool_shared", "jgi_shared", "jgi_exvivo"]
CORI_KNL_QOS = "regular"
CORI_CHARGE_ACCNT = "fungalp"  # cori needs -q, -A
CORI_KNL_CHARGE_ACCNT = "m342"  # for general cori nodes and for KNL
CORI_CONSTRAINT = "haswell"  # [haswell, knl, skylake]
CTR = 0.2  # cloning time rate = life_left / runtime_requested
# ex) if -t=5min, fire sbatch at 4min

CORI_SKYLAKE_CHARGE_ACCNT = "gtrqc"
CORI_SKYLAKE_QOS = "jgi_exvivo"

# ================================
# Lab IT Lawrencium & JGI Cluster
# ================================
LAWRENCIUM_PARTITION = "lr2"
LAWRENCIUM_QOS = "lr_debug"
LAWRENCIUM_CHARGE_ACCNT = "pc_jaws"

JGI_CLOUD_PARTITION = "lr3"
JGI_CLOUD_QOS = "condo_jgicloud"
JGI_CLOUD_ACCNT = "lr_jgicloud"


# =========================
# OLCF Rhea, Summit, Slate
# =========================
RHEA_CHARGE_ACCNT = "bif113"
