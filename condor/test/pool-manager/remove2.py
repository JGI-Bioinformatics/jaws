import os
import time
import subprocess
import shlex
from utils import run_sh_command
import configparser
from datetime import datetime
from operator import itemgetter

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

#
# System vars, dirs, and cmds
#
config = configparser.ConfigParser()
config.read("/global/cfs/cdirs/jaws/condor/pool-manager/jaws_condor_pool_manager.ini")
accnt_name = config["SLURM"]["accnt_name"]
condor_root = config["CONDOR"]["condor_root"]
min_pool_size = config.getint("CONDOR", "min_pool_size")
site_id = config["SLURM"]["site_id"]

def run_scancel(node: str, accnt_name:str, module_load=None):
    cmd = f"""{module_load} squeue --format="%.18i %.40j %R" --me -u {accnt_name} -t R,PD | grep jaws_condor"""
    so, se, ec = run_sh_command(cmd)
    if ec != 0:
        print(f"ERROR: failed to execute command: {cmd}")
        exit(1)
    for l in so.split("\n"):
        if l and len(shlex.split(l)) == 3:
            slurm_jid = int(l.split()[0])
            node_n = l.split()[2]
            if node_n == node:
                cmd = f"{module_load} scancel {slurm_jid}"
                print(cmd)

cmd = "condor_status -autoformat Name  TotalTimeUnclaimedIdle"
so, se, ec = run_sh_command(cmd)
if ec != 0:
    print(f"ERROR: failed to execute command: {cmd}")
    exit(1)
machine_idle_list = []
machine_inuse_list = []
for l in so.split("\n"):
    if l and len(shlex.split(l)) == 2:
        node_name = l.split()[0].split('@')[1]
        idle_time_sec = l.split()[1]
        if l.split()[1] != "undefined":
            machine_idle_list.append([int(idle_time_sec), node_name])
        else:
            if node_name not in machine_inuse_list:
                machine_inuse_list.append(node_name)
machine_idle_list = sorted(machine_idle_list, key=itemgetter(0))
print(machine_idle_list)
print(machine_inuse_list)
for m in machine_idle_list[min_pool_size:]:
    if m[1] not in machine_inuse_list:
        print(m)
        run_scancel(m[1], accnt_name)
        if site_id == "cori": run_scancel(m[1], accnt_name, "module load esslurm &&")