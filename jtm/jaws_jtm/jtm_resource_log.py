#! /usr/bin/env python
#
# Get the task resource consumption log for a given task id
#
# jtm-resource-log returns JSON string if resource log for the task id
#                          None if no resource log
#
import json
import csv
import sys
import os

from jaws_jtm.lib.jtminterface import JtmInterface
from jaws_jtm.lib.common import eprint


def resource_log():
    assert len(sys.argv) == 2, "USAGE: jtm-resource-log <task_id>"
    taskId = int(sys.argv[1])
    resource_log_file = JtmInterface("resource", info_tag=taskId).call(task_id=taskId)
    if resource_log_file == -88:
        eprint("jtm-resource-log: command timeout.")
        sys.exit(-1)
    # print resource_log_file
    # http://www.andymboyle.com/2011/11/02/quick-csv-to-json-parser-in-python/
    if os.path.isfile(resource_log_file):
        f = open(resource_log_file, "rU")
        reader = csv.DictReader(
            f,
            fieldnames=(
                "child_pid",  # 1
                "clone_time_rate",  # 2
                "cpu_load",  # 3
                "end_date",  # 4
                "host_name",  # 5
                "ip_address",  # 6
                "job_time",  # 7
                "life_left",  # 8
                "mem_per_core",  # 9
                "mem_per_node",  # 10
                "num_cores",  # 11
                "num_tasks",  # 12
                "num_workers_on_node",  # 13
                "perc_mem_used",  # 14
                "pool_name",  # 15
                "ret_msg",  # 16
                "rmem_usage",  # 17
                "root_pid",  # 18
                "run_time",  # 19
                "slurm_jobid",  # 20
                "task_id",  # 21
                "vmem_usage",  # 22
                "worker_id",  # 23
                "worker_type",  # 24
                "jtm_host_name",  # 25
                "nwpn",  # 26
            ),
        )
        print(
            """{ "task_id": %d, "resource_log": %s }"""
            % (taskId, json.dumps([row for row in reader]))
        )
    else:
        resource_log_file = None
        print("Resource file, %s, not found." % (resource_log_file))

    sys.exit(0) if resource_log_file is not None else sys.exit(1)
