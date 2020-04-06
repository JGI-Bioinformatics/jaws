import configparser
import os

DEFAULT_CONFIG_FILE = '~/jtm.ini'


# -------------------------------------------------------------------------------------------
# CONSTANTS
# -------------------------------------------------------------------------------------------
class JtmConstants():
    VERSION = "5.6.8"

    # Supported cluster
    COMPUTE_RESOURCES = ["cori",
                         "lawrencium", "jgi_cloud", "jaws_lbl_gov", "jgi_cluster",  "lbl",  # lbl it
                         "rhea", "summit",  # olcf
                         "ssul-dm_dhcp_lbl_gov", "sjsul-lm-2_local", "sjsul-lm.dhcp.lbnl.us",  # testing
                         "summit", "rhea", "slate",  # OLCF
                         "aws"  # future support
                         ]

    # task types
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
              "jtm_host_name": 25,  # hostname on which jtm manager is running ["cori"]
              "nwpn": 26  # num workers per node, NOT USED
              }

    # Note: to see all the qos assigned
    # $ sacctmgr show assoc user=jaws_jtm -p
    # Cluster|Account|User|Partition|Share|Priority|GrpJobs|GrpTRES|GrpSubmit|GrpWall|GrpTRESMins|MaxJobs|MaxTRES|
    # MaxTRESPerNode|MaxSubmit|MaxWall|MaxTRESMins|QOS|Def QOS|GrpTRESRunMins|
    # escori|m342|jaws_jtm||1||||||||bb/datawarp=52828800M|||||debug_hsw,debug_knl,flex,interactive,jupyter,long,low_knl,
    # overrun,premium,regular_0,regular_1,regular_bigmem,resv,resv_shared,shared,xfer|||
    # cori|m342|jaws_jtm||1||||||||bb/datawarp=52828800M|||||debug_hsw,debug_knl,flex,interactive,jupyter,long,low_knl,
    # overrun,premium,regular_0,regular_1,regular_bigmem,resv,resv_shared,shared,xfer|||
    QOS_LIST = ["genepool_special",
                "genepool_shared",
                "genepool",
                "regular",
                "jgi_shared",
                "jgi_exvivo",
                "condo_jgicloud"
                ]
    CONSTRAINT_LIST = ["haswell", "knl", "skylake"]
    DEFAULT_POOL = ["small", "medium", "large", "xlarge"]
    NUM_MANAGER_PROCS = 7
    NUM_WORKER_PROCS = 6


class JtmConfig(object):
    '''
    JTM config utility
    '''
    def __init__(self, config_file=None):
        if not config_file:
            config_file = self.get_config_file()
        self.configparser = self.create_config(config_file)
        self.constants = JtmConstants()

    def get_config_file(self):
        return os.environ.get('JTM_CONFIG_FILE', DEFAULT_CONFIG_FILE)

    def create_config(self, config_file=None):
        config_parser = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
        config_parser.read(config_file or DEFAULT_CONFIG_FILE)
        return config_parser
