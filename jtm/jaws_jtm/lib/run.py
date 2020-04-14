#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)

"""
Created on Apr 3, 2013

Convenience methods to make it easier to run external programs
and other os-related tools

@author: Seung-Jin Sul

"""

from subprocess import Popen, call, PIPE
import os
import glob
import shlex
import re
import time
from threading import Timer
import sys


# -------------------------------------------------------------------------------
def run(*popenargs, **kwargs):
    """
    Run user command using subprocess.call
    :param popenargs: command and options to run
    :param kwargs: additional parameters
    :return:
    """
    kw = {}
    kw.update(kwargs)
    dry_run = kw.pop("dry_run", False)

    if dry_run:
        print(popenargs)
    else:
        # convert something like run("ls -l") into run("ls -l", shell=True)
        if isinstance(popenargs[0], str) and len(shlex.split(popenargs[0])) > 1:
            kw.setdefault("shell", True)

        # > /dev/null 2>&1
        if kw.pop("supressAllOutput", False):
            stdnull = open(os.devnull, "w")  # incompat with close_fds on Windows
            kw.setdefault("stdout", stdnull)
            kw.setdefault("stderr", stdnull)
        else:
            stdnull = None

        return_code = call(*popenargs, **kw)
        if stdnull:
            stdnull.close()
        if return_code != 0:
            eprint("Failed to call run()")
            return 1

    return 0


# -------------------------------------------------------------------------------
# def back_ticks(*popenargs, **kwargs):
#     """
#     Similar to shell backticks, e.g. a = `ls -1` <=> a = backticks(['ls','-1']).
#     If 'dry_run=True' is given as keyword argument, then 'dryRet' keyword must
#     provide a value to return from this function.
#     :param popenargs: command and options to run
#     :param kwargs: additional parameters
#     :return: command result (stdout)
#     """
#     kw = {}
#     kw.update(kwargs)
#     dry_run = kw.pop("dry_run", False)
#     dryRet = kw.pop("dryRet", None)
#
#     if dry_run:
#         print(popenargs)
#         return dryRet
#     else:
#         kw["stdout"] = PIPE
#         p = Popen(*popenargs, **kw)
#         retOut = p.communicate()[0]
#         if p.returncode != 0:
#             eprint("Failed to run back_ticks()")
#             return 1
#
#         return retOut


# -------------------------------------------------------------------------------
def make_dir(path, dry_run=False):
    """
    Create one dir with pathname path or do nothing if it already exists. Same as Linux 'mkdir -p'.
    :param path: path to create
    :param dry_run: dryrun directive
    :return:
    """
    if not dry_run:
        if not os.path.exists(path):
            try:
                original_umask = os.umask(0)
                os.makedirs(path, mode=0o775)
            except Exception:
                pass
            finally:
                os.umask(original_umask)
    else:
        print("make_dir %s" % (path,))

    return 0


# -------------------------------------------------------------------------------
def make_dirs(paths, dry_run=False):
    """
    Create muiltiple dirs with the same semantics as make_dir
    :param paths: path to delete
    :param dry_run: dryrun directive
    :return:
    """
    for path in paths:
        make_dir(path=path, dry_run=dry_run)


# -------------------------------------------------------------------------------
def make_file_path(fileName):
    """
    Assume that the argument is a file name and make all directories that are
    part of it
    :param fileName: create dir to the file
    :return:
    """
    dirName = os.path.dirname(fileName)
    if dirName not in ("", "."):
        make_dir(dirName)


# -------------------------------------------------------------------------------
def rm_dir(path, dry_run=False):
    """
    Remove dir
    :param path: path to delete
    :param dry_run: dryrun directive
    :return:
    """
    # To do: perhaps use shutil.rmtree instead?
    return run(["rm", "-rf", path], dry_run=dry_run)


# -------------------------------------------------------------------------------
def rmf(path, dry_run=False):
    """
    Remove file.
    :param path: path to delete
    :param dry_run: dryrun directive
    :return:
    """
    for f in glob.iglob(path):
        try:
            if os.path.exists(f):
                os.remove(f)
        except OSError:
            pass


# -------------------------------------------------------------------------------
def rmf_many(paths, dry_run=False):
    """
    Remove multiple files.
    :param paths: path to delete
    :param dry_run: dryrun directive
    :return:
    """
    for f in paths:
        try:
            os.remove(f)
        except OSError:
            pass


# -------------------------------------------------------------------------------
def remake_dir(path, dry_run=False):
    """
    Create an empty dir with a given path.
    If path already exists,  it will be removed first.
    :param path: path to delete and create
    :param dry_run:
    :return:
    """
    rm_dir(path, dry_run=dry_run)
    make_dir(path, dry_run=dry_run)


# -------------------------------------------------------------------------------
def chmod(path, mode, opts="", dry_run=False):
    """
    Change mode.
    :param path: path to chmod
    :param mode: the form `[ugoa]*([-+=]([rwxXst]*|[ugo]))+'.
    :param opts: additional chmod options
    :param dry_run: dryrun directive
    :return:
    """
    if isinstance(path, str):
        path = [path]
    else:
        path = list(path)
    run(["chmod"] + opts.split() + [mode] + path, dry_run=dry_run)


# -------------------------------------------------------------------------------
def run_sh_command(cmd, live=True, log=None, run_time=False,
                   show_stdout=True, timeout_sec=0):
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

            # p.wait()

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


# -------------------------------------------------------------------------------
def file_exist(path):
    """
    Check file existence
    """
    return os.path.exists(path)


# -------------------------------------------------------------------------------
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


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
