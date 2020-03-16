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

from jaws_jtm_config import JTM_WORKER_NUM_THREADS
from jaws_jtm.common import eprint, pprint, CalledProcessError, logger
from jaws_jtm.lib.run import back_ticks

g_scale_inv = ((1024.0 * 1024.0, "MB"), (1024.0, "KB"))


def get_cpu_load(pid):
    """
    return cpu usage of process
    :param pid: process id for running "ps"
    :return:
    """
    ps_cmd = "ps h -o pcpu -p %d" % (pid)
    cpu_load = 0

    try:
        ps_stdout_str = back_ticks(ps_cmd, shell=True)
        if sys.platform.lower() == "darwin":
            ps_stdout_str = ps_stdout_str.decode().strip().split("\n")[1]
        cpu_load = ps_stdout_str.strip()
    except CalledProcessError as e:
        eprint("get_cpu_load(): %s" % e)
        raise
    except Exception as ex:
        eprint("get_cpu_load(): %s" % ex)
        raise

    return cpu_load


# -------------------------------------------------------------------------------
def get_runtime(pid):
    """
    get the total runtime in sec for a given pid.
    :param pid: process id for procstat
    :return:
    """
    proc_stat_file = "/proc/%d/stat" % pid
    grep_cmd = 'grep btime /proc/stat | cut -d " " -f 2'
    cat_cmd = 'cat /proc/%d/stat | cut -d " " -f 22' % pid
    proc_run_time = 0

    try:
        boot_time = back_ticks(grep_cmd, shell=True)
        try:
            boot_time = int(boot_time.strip())
        except ValueError:
            boot_time = 0

        if os.path.exists(proc_stat_file):
            msec_since_boot = back_ticks(cat_cmd, shell=True)
        else:
            return proc_run_time

        try:  # To handle ValueError: invalid literal for int() with base 10: ''
            msec_since_boot = int(msec_since_boot.strip())
            sec_since_boot = msec_since_boot / 100
        except ValueError:
            sec_since_boot = 0

        proc_start_time = boot_time + sec_since_boot
        now = time.time()
        proc_run_time = int(now - proc_start_time)  # unit in seconds will be enough
    except CalledProcessError as msg:
        logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))

    return proc_run_time


# -------------------------------------------------------------------------------
# def get_memory_usage(pid):
#    """
#    return memory usage in Mb.
#
#    :param pid:
#    """
#    return _VmB('VmSize:', pid)
# -------------------------------------------------------------------------------
def _VmB(VmKey, pid):
    """
    get various mem usage properties of process with id pid in MB
    :param VmKey:
    :param pid:
    :return:
    """
    proc_status_loc = "/proc/%d/status" % pid
    unit_scale = {"kB": 1.0 / 1024.0, "mB": 1.0, "KB": 1.0 / 1024.0, "MB": 1.0}

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
        # logger.error("VmKey, {}, not found in {}".format(VmKey, v))
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
def toScale(x):
    """
    convert scale: convert 'x' to a string with B/KB/MB units
    :param x: integer or float
    :return:

    >>> toScale("no int")
    Traceback (most recent call last):
        ...
    TypeError: unsupported operand type(s) for /: 'str' and 'float'
    >>> toScale(100.0)
    '0.098B'
    >>> toScale(1024)
    '1.000KB'
    >>> toScale(1024*1024)
    '1.000MB'
    """
    for sc in g_scale_inv:
        y = x / sc[0]
        if y >= 1:
            return "%.3f%s" % (y, sc[1])
    return "%.3f%s" % (y, "B")


# -------------------------------------------------------------------------------
def get_virtual_memory_usage(pid, since=0.0, as_str=True):
    """
    Return memory usage in bytes or as formatted string.
    :param pid:
    :param since:
    :param as_str:
    :return:
    """
    b = _VmB("VmSize:", pid) - since
    if as_str:
        return "VirtMem: " + toScale(b)
    else:
        return b


# -------------------------------------------------------------------------------
def get_resident_memory_usage(pid, since=0.0, as_str=True):
    """
    Return resident memory usage in bytes.
    :param pid:
    :param since:
    :param as_str:
    :return:

    >>> get_resident_memory_usage(0)
    'ResMem: 0.000B'
    """
    b = _VmB("VmRSS:", pid) - since
    if as_str:
        return "ResMem: " + toScale(b)
    else:
        return b


# -------------------------------------------------------------------------------
def get_stacksize(pid, since=0.0, as_str=True):
    """
    Return stack size in bytes.
    :param pid:
    :param since:
    :param as_str:
    :return:
    """
    b = _VmB("VmStk:", pid) - since
    if as_str:
        return "StackMem: " + toScale(b)
    else:
        return b


# -------------------------------------------------------------------------------
def get_pid_tree(pid):
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
        import psutil

        try:
            child_pid_list.extend(
                [p.pid for p in psutil.Process(pid).children(recursive=True)]
            )
        except psutil.NoSuchProcess:
            logger.exception(
                "Failed to call psutil.Process(). Process id is not exist."
            )
        except Exception as psutil_error:
            logger.exception(psutil_error)
            raise
    else:
        cmd = "ps -o pid --ppid %d --noheaders" % pid
        child_pid_list.append(pid)
        try:
            # ps_stdout_str = back_ticks(cmd, shell=True)
            ps_stdout_str = subprocess.Popen(
                [cmd], shell=True, stdout=subprocess.PIPE
            ).communicate()[0]
            child_pid_list.extend(
                [int(pidStr) for pidStr in ps_stdout_str.split("\n")[:-1]]
            )
        except CalledProcessError as msg:
            logger.exception(
                "Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode)
            )
            child_pid_list = []
        except Exception as ps_error:
            logger.exception(ps_error)
            # raise

    return child_pid_list


# -------------------------------------------------------------------------------
def get_total_mem_usage_per_node():
    """
    get % mem used per node
    :return:
    """
    mem_perc = 0.0
    if sys.platform.lower() == "darwin":
        cmd = "top -l 1 | head -n 10 | grep PhysMem | sed 's/,//g'"
        out = back_ticks(cmd, shell=True)
        # logger.debug(out)
        p = r"(\d+)"
        f = re.findall(p, out.decode())
        mem_perc = 100.0 - float(f[2]) / (float(f[0]) * 1024) * 100.0
    else:  # linux
        cmd = "free"
        try:
            freeOut = back_ticks(cmd, shell=True)
        except CalledProcessError as msg:
            logger.exception(
                "Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode)
            )
            return -1
        mem_perc = (
            float(freeOut.split("\n")[1].split()[2])
            / float(freeOut.split("\n")[1].split()[1])
            * 100.0
        )

    return mem_perc


# -------------------------------------------------------------------------------
def get_num_workers_on_node():
    """
    get total number of workers on a given node
    :return:

    >>> get_num_workers_on_node()
    0
    """
    cmd = "ps ax | grep -v grep | grep jtm-worker | wc -l"
    # ps_stdout_str = 0

    try:
        ps_stdout_str = back_ticks(cmd, shell=True)
    except CalledProcessError as msg:
        logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
        return -1

    return (
        int(ps_stdout_str) / JTM_WORKER_NUM_THREADS
    )  # we have three processes per a worker


# -------------------------------------------------------------------------------
def darwin_free():
    """
    mac os version free()
    :return:
    """
    # Get process info
    ps = (
        subprocess.Popen(["ps", "-caxm", "-orss,comm"], stdout=subprocess.PIPE)
        .communicate()[0]
        .decode()
    )
    vm = subprocess.Popen(["vm_stat"], stdout=subprocess.PIPE).communicate()[0].decode()

    # Iterate processes
    ps_stdout_split_list = ps.split("\n")
    sep = re.compile(r'[\s]+')
    rss_total = 0  # kB
    for row in range(1, len(ps_stdout_split_list)):
        row_text = ps_stdout_split_list[row].strip()
        row_element = sep.split(row_text)
        try:
            rss = float(row_element[0]) * 1024
        except:
            rss = 0  # ignore...
        rss_total += rss

    # Process vm_stat
    vm_lines_list = vm.split("\n")
    sep = re.compile(r':[\s]+')
    vm_stats_dict = {}
    for row in range(1, len(vm_lines_list) - 2):
        row_text = vm_lines_list[row].strip()
        row_element = sep.split(row_text)
        vm_stats_dict[(row_element[0])] = int(row_element[1].strip("\.")) * 4096

    pprint.pprint(vm_stats_dict)
    print(
        ("Wired Memory:\t\t%d MB" % (vm_stats_dict["Pages wired down"] / 1024 / 1024))
    )
    print(("Active Memory:\t\t%d MB" % (vm_stats_dict["Pages active"] / 1024 / 1024)))
    print(("Inactive Memory:\t%d MB" % (vm_stats_dict["Pages inactive"] / 1024 / 1024)))
    print(("Free Memory:\t\t%d MB" % (vm_stats_dict["Pages free"] / 1024 / 1024)))
    print(("Real Mem Total (ps):\t%.3f MB" % (rss_total / 1024 / 1024)))
