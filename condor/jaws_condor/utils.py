#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)

"""
Created on Apr 3, 2013

Convenience methods to make it easier to run external programs
and other os-related tools

@author: Seung-Jin Sul

"""

from subprocess import Popen, PIPE
import re
import time
from threading import Timer


def run_sh_command(
    cmd, live=True, log=None, run_time=False, show_stdout=True, timeout_sec=0
):
    """
    Run a command, catch stdout and stderr and exit_code
    :param cmd:
    :param live: live (boolean, default True - don't run the command but pretend we did)
    :param log:
    :param run_time: flag to print elapsed time as log
    :param show_stdout: flag to show stdout or not
    :param timeout_sec: timeout to terminate the specified command
    :return:

    >>> run_sh_command("ls -al", live=False)
    ("Not live: cmd = 'ls -al'", None, 0)
    """
    std_out = None
    std_err = None
    exit_code = None
    start = 0
    end = 0
    elapsed_sec = 0

    if cmd:
        if not live:
            std_out = "Not live: cmd = '%s'" % (cmd)
            exit_code = 0
        else:
            if log and show_stdout:
                log.info("cmd: %s" % (cmd))

            # ---------
            # OLD
            if run_time:
                start = time.time()

            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)

            # ref) http://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
            if timeout_sec > 0:
                timer = Timer(timeout_sec, p.kill)

            try:
                if timeout_sec > 0:
                    timer.start()
                std_out, std_err = p.communicate()
                if type(std_out) is bytes:
                    std_out = std_out.decode()
                if type(std_err) is bytes:
                    std_err = std_err.decode()
                if show_stdout:
                    print(std_out)  # for printing slurm job id
                exit_code = p.returncode
            finally:
                if timeout_sec > 0:
                    timer.cancel()
                    exit_code = 143

            if run_time:
                end = time.time()
                elapsed_sec = end - start
                if log:
                    log.info("*************************************")
                    if cmd.split(" ")[0].split("/")[-1]:
                        log.info(cmd.split(" ")[0].split("/")[-1])
                    log.info("Command took " + str(elapsed_sec) + " sec.")

                    log.info("*************************************")

            if log and show_stdout:
                log.info(
                    "Return values: exit_code="
                    + str(exit_code)
                    + ", std_out="
                    + str(std_out)
                    + ", std_err="
                    + str(std_err)
                )

            if exit_code != 0:
                if log:
                    log.warn("- The exit code has non-zero value.")

    else:
        if log:
            log.error("- No command to run.")
            return None, None, -1

    if type(std_out) is bytes:
        std_out = std_out.decode()
    if type(std_err) is bytes:
        std_err = std_err.decode()

    return std_out, std_err, exit_code


def pad_string_path(my_string, pad_length=8, depth=None):
    """
    Pads a string with 0's and splits into groups of 2
    e.g. pad_stringPath(44) returns "00/00/00/44"

    @param my_string: string to pad
    @return: padded string
    """
    my_string = str(my_string)
    pad_length = int(pad_length)

    if pad_length > 8 or pad_length <= 0:
        pad_length = 8

    # left-pad with 0's
    my_string = my_string.zfill(pad_length)

    # use re.findall function to split into pairs of strings
    string_list = re.findall("..", my_string)

    # create ss/ss/ss/ss
    if not depth:
        pad_string = "/".join(string_list)
    else:
        pad_string = "/".join(string_list[:depth])

    pad_string = pad_string + "/"

    return pad_string


def extract_cromwell_id(task: str) -> str:
    """
    Extract Cromwell run id from a Cromwell task string

    ex)
    /bin/bash /..../cromwell-executions/x/74a668ad-958e-437e-a941-6e5e23f8716d/call-y/shard-01/execution/script
    ==> extract 74a668ad-958e-437e-a941-6e5e23f8716d

    :param task: command string
    :return: Cromwell run id in UUID format
    """
    cromwell_id = None
    # NOTE: here UUID format spec is assumed to comply with "8-4-4-4-12" format
    try:
        regex = re.compile(r"cromwell-executions\/[^\/]+\/([^\/]+)", re.I)
        match = regex.search(task)
    except Exception:
        raise

    if match:
        cromwell_id = match.group(1)

    return cromwell_id


def run_slurm_cmd(s_cmd: str) -> int:
    """
    To run sbatch and squeue and parse `wc -l`
    """
    so, se, ec = run_sh_command(s_cmd, show_stdout=False)
    print(s_cmd)
    if ec != 0:
        print(f"ERROR: failed to execute slurm command: {s_cmd}\n{se}")
        exit(1)
    try:
        return int(so.rstrip())
    except TypeError as te:
        print(f"Unexpected output from slurm command: {s_cmd}\n{te}\n{se}")
        return -1


def mem_unit_to_g(unit: str, mem: float) -> float:
    if unit.upper() not in ("GB", "G"):
        if unit.upper() in ("TB", "T"):
            mem *= 1024
        elif unit.upper() in ("MB", "M"):
            mem /= 1024
        elif unit.upper() in ("KB", "K"):
            mem /= 1024 * 1024
        else:
            return -1
    return mem
