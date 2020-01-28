#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Copyright 2018-2019 (C) DOE JGI, LBL
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

from Common import *
from Utils.Run import *

g_scale_inv = ( (1024.*1024., "MB"), (1024., "KB") )

#-------------------------------------------------------------------------------
def get_cpu_load(pid):
    """
    return cpu usage of process
    :param pid: process id for running "ps"
    :return:
    """
    psCmd = "ps h -o pcpu -p %d" % (pid)
    cpuLoad = 0

    try:
        psOut = backTicks(psCmd, shell=True)
        if sys.platform.lower() == "darwin":
            psOut = psOut.strip().split('\n')[1]

        cpuLoad = psOut.strip()
    except CalledProcessError:
        #logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
        pass

    return cpuLoad


#-------------------------------------------------------------------------------
def get_runtime(pid):
    """
    get the total runtime in sec for a given pid.
    :param pid: process id for procstat
    :return:
    """
    procStatFile = '/proc/%d/stat' % pid
    grepCmd = 'grep btime /proc/stat | cut -d " " -f 2'
    catCmd = 'cat /proc/%d/stat | cut -d " " -f 22' % pid
    procRunTime = 0

    try:
        bootTime = backTicks(grepCmd, shell=True)
        try:
            bootTime = int(bootTime.strip())
        except ValueError:
            bootTime = 0

        if os.path.exists(procStatFile):
            msecSinceBoot = backTicks(catCmd, shell=True)
        else:
            return procRunTime

        try:  # To handle ValueError: invalid literal for int() with base 10: ''
            msecSinceBoot = int(msecSinceBoot.strip())
            secSinceBoot = msecSinceBoot / 100
        except ValueError:
            # msecSinceBoot = 0
            secSinceBoot = 0

        pStartTime = bootTime + secSinceBoot
        now = time.time()
        procRunTime = int(now - pStartTime)  # unit in seconds will be enough
    except CalledProcessError as msg:
        logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
        pass

    return procRunTime


##-------------------------------------------------------------------------------
#def get_memory_usage(pid):
#    """
#    return memory usage in Mb.
#
#    :param pid:
#    """
#    return _VmB('VmSize:', pid)


#-------------------------------------------------------------------------------
def _VmB(VmKey, pid):
    """
    get various mem usage properties of process with id pid in MB
    :param VmKey:
    :param pid:
    :return:
    """
    procStatus = '/proc/%d/status' % pid
    unitScale = {'kB': 1.0/1024.0, 'mB': 1.0,
                 'KB': 1.0/1024.0, 'MB': 1.0}

    ## get pseudo file /proc/<pid>/status
    try:
        if os.path.exists(procStatus):
            t = open(procStatus)
            v = t.read()
            t.close()
        else:
            return 0.0
    except OSError:
        logger.exception("Failed to open /proc files.")
        return 0.0  # non-Linux?

    ## get VmKey line e.g. 'VmRSS: 9999 kB\n ...'
    try:
        i = v.index(VmKey)
        v = v[i:].split(None, 3)  # by whitespace
    except (ValueError, UnboundLocalError):
        # logger.error("VmKey, {}, not found in {}".format(VmKey, v))
        pass

    if len(v) < 3:
        return 0.0  # invalid format?

    ## convert Vm value to bytes
    try:
        byteVal = float(v[1]) * unitScale[v[2]]
    except ValueError:
        byteVal = 0.0


    return byteVal


#-------------------------------------------------------------------------------
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
        y = x/sc[0]
        if y >= 1:
            return "%.3f%s" % (y, sc[1])
    return "%.3f%s" % (y, "B")


#-------------------------------------------------------------------------------
def get_virtual_memory_usage(pid, since=0.0, asStr=True):
    """
    Return memory usage in bytes or as formatted string.
    :param pid:
    :param since:
    :param asStr:
    :return:
    """
    b = _VmB('VmSize:', pid) - since
    if asStr:
        return "VirtMem: " + toScale(b)
    else:
        return b


#-------------------------------------------------------------------------------
def get_resident_memory_usage(pid, since=0.0, asStr=True):
    """
    Return resident memory usage in bytes.
    :param pid:
    :param since:
    :param asStr:
    :return:

    >>> get_resident_memory_usage(0)
    'ResMem: 0.000B'
    """
    b = _VmB('VmRSS:', pid) - since
    if asStr:
        return "ResMem: " + toScale(b)
    else:
        return b


#-------------------------------------------------------------------------------
def get_stacksize(pid, since=0.0, asStr=True):
    """
    Return stack size in bytes.
    :param pid:
    :param since:
    :param asStr:
    :return:
    """
    b = _VmB('VmStk:', pid) - since
    if asStr:
        return "StackMem: " + toScale(b)
    else:
        return b


#-------------------------------------------------------------------------------
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
    cpids = []
    if sys.platform.lower() == "darwin":  # if mac os
        ## Todo
        import psutil
        try:
            cpids.extend([p.pid for p in psutil.Process(pid).children(recursive=True)])
        except psutil.NoSuchProcess:
            logger.exception("Failed to call psutil.Process().")
            pass
        except Exception as psutil_error:
            logger.exception(psutil_error)
            raise
            # pass
    else:
        cmd = "ps -o pid --ppid %d --noheaders" % pid
        cpids.append(pid)
        try:
            # psOut = backTicks(cmd, shell=True)
            psOut = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE).communicate()[0]
            cpids.extend([int(pidStr) for pidStr in psOut.split("\n")[:-1]])
        except CalledProcessError as msg:
            logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
            cpids = []
            pass
        except Exception as ps_error:
            logger.exception(ps_error)
            # raise
            pass

    return cpids


#-------------------------------------------------------------------------------
def get_total_mem_usage_per_node():
    """
    get % mem used per node
    :return:
    """
    memPerc = 0.0
    if sys.platform.lower() == "darwin":
        cmd = "top -l 1 | head -n 10 | grep PhysMem | sed 's/,//g'"
        out = backTicks(cmd, shell=True)
        # logger.debug(out)
        p = r"(\d+)"
        f = re.findall(p, out)
        memPerc = 100.0 - float(f[2]) / (float(f[0]) * 1024) * 100.0
    else:  # linux
        cmd = "free"
        try:
            freeOut = backTicks(cmd, shell=True)
        except CalledProcessError as msg:
            logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
            return -1
            memPerc = float(freeOut.split('\n')[1].split()[2]) / \
                      float(freeOut.split('\n')[1].split()[1]) * 100.0

    return memPerc


#-------------------------------------------------------------------------------
def get_num_workers_on_node():
    """
    get total number of workers on a given node
    :return:

    >>> get_num_workers_on_node()
    0
    """
    cmd = "ps ax | grep -v grep | grep jtm-worker | wc -l"
    # psOut = 0

    try:
        psOut = backTicks(cmd, shell=True)
    except CalledProcessError as msg:
        logger.exception("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
        return -1

    return int(psOut) / 3 ## we have three processes per a worker


#-------------------------------------------------------------------------------
def darwin_free():
    """
    mac os version free()
    :return:
    """
    # Get process info
    ps = subprocess.Popen(['ps', '-caxm', '-orss,comm'], stdout=subprocess.PIPE).communicate()[0].decode()
    vm = subprocess.Popen(['vm_stat'], stdout=subprocess.PIPE).communicate()[0].decode()

    # Iterate processes
    processLines = ps.split('\n')
    sep = re.compile('[\s]+')
    rssTotal = 0  # kB
    for row in range(1, len(processLines)):
        rowText = processLines[row].strip()
        rowElements = sep.split(rowText)
        try:
            rss = float(rowElements[0]) * 1024
        except:
            rss = 0  # ignore...
        rssTotal += rss

    # Process vm_stat
    vmLines = vm.split('\n')
    sep = re.compile(':[\s]+')
    vmStats = {}
    for row in range(1, len(vmLines) - 2):
        rowText = vmLines[row].strip()
        rowElements = sep.split(rowText)
        vmStats[(rowElements[0])] = int(rowElements[1].strip('\.')) * 4096

    pprint.pprint(vmStats)
    print('Wired Memory:\t\t%d MB' % (vmStats["Pages wired down"] / 1024 / 1024))
    print('Active Memory:\t\t%d MB' % (vmStats["Pages active"] / 1024 / 1024))
    print('Inactive Memory:\t%d MB' % (vmStats["Pages inactive"] / 1024 / 1024))
    print('Free Memory:\t\t%d MB' % (vmStats["Pages free"] / 1024 / 1024))
    print('Real Mem Total (ps):\t%.3f MB' % (rssTotal / 1024 / 1024))


#-------------------------------------------------------------------------------
if __name__ == "__main__":
    # doctest: python ResourceUsage.py -v
    import doctest
    doctest.testmod()