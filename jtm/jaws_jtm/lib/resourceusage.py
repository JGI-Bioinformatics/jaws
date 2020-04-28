#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)

"""
Get resource usage for the process.
Adapted from http://code.activestate.com/recipes/286222/
"""

import os
import sys
import time
import subprocess
import re
import psutil
import logging

from jaws_jtm.lib.run import run_sh_command
from jaws_jtm.common import logger

SCALE_INV = ((1024.*1024., "MB"), (1024., "KB"))
logger = logging.getLogger(__package__)


# -------------------------------------------------------------------------------
def get_cpu_load(pid: int) -> float:
    """
    return cpu usage of process
    :param pid: process id for running "ps"
    :return:
    """
    ps_cmd = "ps h -o pcpu -p %d" % (pid)
    cpu_load = 0.0

    try:
        ps_stdout_str, _, _ = run_sh_command(ps_cmd, log=logger, show_stdout=False)
        if sys.platform.lower() == "darwin":
            ps_stdout_str = ps_stdout_str.strip().split('\n')[1]
        cpu_load = float(ps_stdout_str.strip())
    except IndexError as e:
        eprint(e)
        cpu_load = 0.0
    except ValueError as e:
        eprint(e)
        cpu_load = 0.0
    except subprocess.CalledProcessError as e:
        eprint("get_cpu_load(): %s" % e)
        raise
    except Exception as ex:
        eprint("get_cpu_load(): %s" % ex)
        raise

    return cpu_load


# -------------------------------------------------------------------------------
def get_runtime(pid: int) -> int:
    """
    get the total runtime in sec for a given pid.
    :param pid: process id for procstat
    :return:
    """
    proc_stat_file = '/proc/%d/stat' % pid
    grep_cmd = 'grep btime /proc/stat | cut -d " " -f 2'
    cat_cmd = 'cat /proc/%d/stat | cut -d " " -f 22' % pid
    proc_run_time = 0

    try:
        boot_time, _, _ = run_sh_command(grep_cmd, log=logger, show_stdout=False)
    except subprocess.CalledProcessError as msg:
        logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))

    try:
        boot_time = int(boot_time.strip())
    except ValueError:
        boot_time = 0

    try:
        if os.path.exists(proc_stat_file):
            msec_since_boot, _, _ = run_sh_command(cat_cmd, log=logger, show_stdout=False)
        else:
            return 0
    except subprocess.CalledProcessError as msg:
        logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))

    try:  # To handle ValueError: invalid literal for int() with base 10: ''
        msec_since_boot = int(msec_since_boot.strip())
        sec_since_boot = msec_since_boot / 100
    except ValueError:
        sec_since_boot = 0

    proc_start_time = boot_time + sec_since_boot
    now = time.time()
    proc_run_time = int(now - proc_start_time)  # unit in seconds will be enough

    return proc_run_time


# -------------------------------------------------------------------------------
def _VmB(VmKey: str, pid: int) -> int:
    """
    get various mem usage properties of process with id pid in MB
    :param VmKey:
    :param pid:
    :return:
    """
    proc_status_loc = '/proc/%d/status' % pid
    unit_scale = {'kB': 1.0/1024.0, 'mB': 1.0,
                  'KB': 1.0/1024.0, 'MB': 1.0}

    # get pseudo file /proc/<pid>/status
    try:
        if os.path.exists(proc_status_loc):
            t = open(proc_status_loc)
            v = t.read()
            t.close()
        else:
            return 0.0
    except OSError:
        logger.exception("Failed to open /proc files.")
        return 0.0  # non-Linux?

    # get VmKey line e.g. 'VmRSS: 9999 kB\n ...'
    try:
        i = v.index(VmKey)
        v = v[i:].split(None, 3)  # by whitespace
    except (ValueError, UnboundLocalError):
        pass

    if len(v) < 3:
        return 0.0  # invalid format?

    # convert Vm value to bytes
    try:
        byte_value = float(v[1]) * unit_scale[v[2]]
    except ValueError:
        byte_value = 0.0

    return byte_value


# -------------------------------------------------------------------------------
def convert_scale(x: int) -> str:
    """
    convert scale: convert 'x' to a string with B/KB/MB units
    :param x: integer or float
    :return:

    >>> convert_scale("no int")
    Traceback (most recent call last):
        ...
    TypeError: unsupported operand type(s) for /: 'str' and 'float'
    >>> convert_scale(100.0)
    '0.098B'
    >>> convert_scale(1024)
    '1.000KB'
    >>> convert_scale(1024*1024)
    '1.000MB'
    """
    for sc in SCALE_INV:
        y = x/sc[0]
        if y >= 1:
            return "%.3f%s" % (y, sc[1])
    return "%.3f%s" % (y, "B")


# -------------------------------------------------------------------------------
def get_virtual_memory_usage(pid: int, since=0.0, as_str=True):
    """
    Return memory usage in bytes or as formatted string.
    :param pid:
    :param since:
    :param as_str:
    :return:
    """
    if sys.platform.lower() == "darwin":
        # NEW
        process = psutil.Process(pid)
        return process.memory_info().vms / 1024.0 / 1024.0  # in MB
    else:
        b = _VmB('VmSize:', pid) - since
        if as_str:
            return "VirtMem: " + convert_scale(b)
        else:
            return b


# -------------------------------------------------------------------------------
def get_resident_memory_usage(pid: int, since=0.0, as_str=True):
    """
    Return resident memory usage in bytes.
    :param pid:
    :param since:
    :param as_str:
    :return:

    >>> get_resident_memory_usage(0)
    'ResMem: 0.000B'
    """
    if sys.platform.lower() == "darwin":
        # NEW
        process = psutil.Process(pid)
        # print("rss: %d" % process.memory_info().rss)
        return process.memory_info().rss / 1024.0 / 1024.0  # in MB
    else:
        b = _VmB('VmRSS:', pid) - since
        if as_str:
            return "ResMem: " + convert_scale(b)
        else:
            return b


# -------------------------------------------------------------------------------
def get_stacksize(pid: int, since=0.0, as_str=True):
    """
    Return stack size in bytes.
    :param pid:
    :param since:
    :param as_str:
    :return:
    """
    b = _VmB('VmStk:', pid) - since
    if as_str:
        return "StackMem: " + convert_scale(b)
    else:
        return b


# -------------------------------------------------------------------------------
def get_pid_tree(pid: int) -> list:
    """
    get the process id tree from parentPid = root to leaf processes
    :param pid:
    :return:

    >>> get_pid_tree(22.12)
    Traceback (most recent call last):
        ...
    TypeError: pid must be an integer (got 22.12)
    >>> get_pid_tree("test")
    Traceback (most recent call last):
        ...
    TypeError: pid must be an integer (got 'test')
    """
    child_pid_list = []
    if sys.platform.lower() == "darwin":  # if mac os
        # Todo
        try:
            child_pid_list.extend([p.pid for p in psutil.Process(pid).children(recursive=True)])
        except (psutil.NoSuchProcess, ProcessLookupError):
            logger.warning("Failed to call psutil.Process(). Process id is not exist.")
        except Exception as psutil_error:
            logger.warning(psutil_error)

    else:
        cmd = "ps -o pid --ppid %d --noheaders" % pid
        child_pid_list.append(pid)
        try:
            ps_stdout_str = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE).communicate()[0]

            if type(ps_stdout_str) is bytes:
                ps_stdout_str = ps_stdout_str.decode()

            child_pid_list.extend([int(pidStr) for pidStr in ps_stdout_str.split("\n")[:-1]])
        except subprocess.CalledProcessError as msg:
            logger.warning("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
            child_pid_list = []
        except Exception as ps_error:
            logger.warning(ps_error)

    return child_pid_list


# -------------------------------------------------------------------------------
def get_total_mem_usage_per_node() -> float:
    """
    get % mem used per node
    :return:
    """
    mem_perc = 0.0
    if sys.platform.lower() == "darwin":
        cmd = "top -l 1 | head -n 10 | grep PhysMem | sed 's/,//g'"
        try:
            top_out, _, _ = run_sh_command(cmd, log=logger, show_stdout=False)
        except subprocess.CalledProcessError as msg:
            logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
            return -1
        p = r"(\d+)"
        f = re.findall(p, top_out)
        try:
            mem_perc = 100.0 - float(f[2]) / (float(f[0]) * 1024) * 100.0
        except ZeroDivisionError:
            mem_perc = 0

    else:  # linux
        cmd = "free"
        try:
            free_output, _, _ = run_sh_command(cmd, log=logger, show_stdout=False)
        except subprocess.CalledProcessError as msg:
            logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
            return -1
        try:
            mem_perc = float(free_output.split('\n')[1].split()[2]) / \
                       float(free_output.split('\n')[1].split()[1]) * 100.0
        except ZeroDivisionError:
            mem_perc = 0

    return mem_perc


# -------------------------------------------------------------------------------
def get_num_workers_on_node(config=None) -> int:
    """
    Get total number of workers on a given node
    :return:

    >>> get_num_workers_on_node()
    0
    """
    cmd = "ps ax | grep -v grep | grep jtm-worker | wc -l"
    ps_stdout_str = None

    try:
        ps_stdout_str, _, _ = run_sh_command(cmd, log=logger, show_stdout=False)
    except subprocess.CalledProcessError as msg:
        logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
        return -1

    NUM_WORKER_PROCS = 6
    if config:
        NUM_WORKER_PROCS = config.constants.NUM_WORKER_PROCS

    return int(ps_stdout_str) / NUM_WORKER_PROCS


# -------------------------------------------------------------------------------
def darwin_free_mem() -> int:
    """
    mac os version free()
    :return:
    """
    # Get process info
    ps = subprocess.Popen(['ps', '-caxm', '-orss,comm'],
                          stdout=subprocess.PIPE).communicate()[0].decode()
    vm = subprocess.Popen(['vm_stat'],
                          stdout=subprocess.PIPE).communicate()[0].decode()

    # Iterate processes
    ps_stdout_split_list = ps.split('\n')
    sep = re.compile(r'[\s]+')
    rss_total = 0  # kB
    for row in range(1, len(ps_stdout_split_list)):
        row_text = ps_stdout_split_list[row].strip()
        row_element = sep.split(row_text)
        try:
            rss = float(row_element[0]) * 1024
        except Exception:
            rss = 0  # ignore...
        rss_total += rss

    # Process vm_stat
    vm_lines_list = vm.split('\n')
    sep = re.compile(r':[\s]+')
    vm_stats_dict = {}
    for row in range(1, len(vm_lines_list) - 2):
        row_text = vm_lines_list[row].strip()
        row_element = sep.split(row_text)
        # print(row_element)
        vm_stats_dict[(row_element[0])] = int(row_element[1].strip('.')) * 4096

    # byte unit
    return vm_stats_dict["Pages free"]


# -------------------------------------------------------------------------------
def mem_rss_usage_pid_darwin() -> int:
    """
    Get memory size for rss
    :return:
    """
    process = psutil.Process(os.getpid())
    return process.memory_info().rss  # in bytes


# -------------------------------------------------------------------------------
def mem_vms_usage_pid_darwin() -> int:
    """
    Get memory size for vms
    :return:
    """
    process = psutil.Process(os.getpid())
    return process.memory_info().vms  # in bytes


# -------------------------------------------------------------------------------
def get_free_memory() -> int:
    """
    Get free memory size for Linux
    :return:
    """
    free_mem_bytes = 0
    if sys.platform.lower() == "darwin":  # if mac os
        free_mem_bytes = darwin_free_mem()
    else:
        with open('/proc/meminfo', 'r') as mem:
            for line in mem:
                sline = line.split()
                if str(sline[0]) == 'MemAvailable:':
                    free_mem_bytes = int(sline[1]) * 1024.0  # Bytes
                    break

    # byte unit
    return free_mem_bytes
