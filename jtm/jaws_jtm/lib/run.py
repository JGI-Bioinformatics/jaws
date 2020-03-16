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
import time
from time import Timer
import re
from jasws_jtm.common import eprint


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
    dryRun = kw.pop("dryRun", False)

    if dryRun:
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

        returncode = call(*popenargs, **kw)
        if stdnull:
            stdnull.close()
        if returncode != 0:
            eprint("Failed to call run()")
            return 1

    return 0


# -------------------------------------------------------------------------------
def back_ticks(*popenargs, **kwargs):
    """
    Similar to shell backticks, e.g. a = `ls -1` <=> a = backticks(['ls','-1']).
    If 'dryRun=True' is given as keyword argument, then 'dryRet' keyword must
    provide a value to return from this function.
    :param popenargs: command and options to run
    :param kwargs: additional parameters
    :return: command result (stdout)
    """
    kw = {}
    kw.update(kwargs)
    dryRun = kw.pop("dryRun", False)
    dryRet = kw.pop("dryRet", None)

    if dryRun:
        print(popenargs)
        return dryRet
    else:
        kw["stdout"] = PIPE
        p = Popen(*popenargs, **kw)
        retOut = p.communicate()[0]
        if p.returncode != 0:
            eprint("Failed to run back_ticks()")
            return 1

        return retOut


# -------------------------------------------------------------------------------
def make_dir(path, dryRun=False):
    """
    Create one dir with pathname path or do nothing if it already exists. Same as Linux 'mkdir -p'.
    :param path: path to create
    :param dryRun: dryrun directive
    :return:
    """
    if not dryRun:
        if not os.path.exists(path):
            try:
                original_umask = os.umask(0)
                os.makedirs(path, mode=0o775)
            except:
                pass
            finally:
                os.umask(original_umask)
    else:
        print(("make_dir %s" % (path,)))

    return 0


# -------------------------------------------------------------------------------
def make_dirs(paths, dryRun=False):
    """
    Create muiltiple dirs with the same semantics as make_dir
    :param paths: path to delete
    :param dryRun: dryrun directive
    :return:
    """
    for path in paths:
        make_dir(path=path, dryRun=dryRun)


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
def rm_dir(path, dryRun=False):
    """
    Remove dir
    :param path: path to delete
    :param dryRun: dryrun directive
    :return:
    """
    # To do: perhaps use shutil.rmtree instead?
    return run(["rm", "-rf", path], dryRun=dryRun)


# make alias
rmrf = rm_dir


# -------------------------------------------------------------------------------
def rmf(path, dryRun=False):
    """
    Remove file.
    :param path: path to delete
    :param dryRun: dryrun directive
    :return:
    """
    for f in glob.iglob(path):
        try:
            if os.path.exists(f):
                os.remove(f)
        except OSError:
            pass


# -------------------------------------------------------------------------------
def rmf_many(paths, dryRun=False):
    """
    Remove multiple files.
    :param paths: path to delete
    :param dryRun: dryrun directive
    :return:
    """
    for f in paths:
        try:
            os.remove(f)
        except OSError:
            pass


# -------------------------------------------------------------------------------
def remake_dir(path, dryRun=False):
    """
    Create an empty dir with a given path.
    If path already exists,  it will be removed first.
    :param path: path to delete and create
    :param dryRun:
    :return:
    """
    rmrf(path, dryRun=dryRun)
    make_dir(path, dryRun=dryRun)


# -------------------------------------------------------------------------------
def chmod(path, mode, opts="", dryRun=False):
    """
    Change mode.
    :param path: path to chmod
    :param mode: the form `[ugoa]*([-+=]([rwxXst]*|[ugo]))+'.
    :param opts: additional chmod options
    :param dryRun: dryrun directive
    :return:
    """
    if isinstance(path, str):
        path = [path]
    else:
        path = list(path)
    run(["chmod"] + opts.split() + [mode] + path, dryRun=dryRun)


# -------------------------------------------------------------------------------
def run_sh_command(
    cmd, live=False, log=None, runTime=False, stdoutPrint=True, timeoutSec=0
):
    """
    Run a command, catch stdout and stderr and exitCode
    :param cmd:
    :param live: live (boolean, default false - don't run the command but pretend we did)
    :param log:
    :param runTime: flag to print elapsed time as log
    :param stdoutPrint: flag to show stdout or not
    :param timeoutSec: timeout to terminate the specified command
    :return:

    >>> run_sh_command("ls -al", live=False)
    ("Not live: cmd = 'ls -al'", None, 0)
    """
    stdOut = None
    stdErr = None
    exitCode = None
    start = 0
    end = 0
    elapsedSec = 0

    if cmd:
        if not live:
            stdOut = "Not live: cmd = '%s'" % (cmd)
            exitCode = 0
        else:
            if log and stdoutPrint:
                log.info("cmd: %s" % (cmd))

            # ---------
            # OLD
            if runTime:
                start = time.time()

            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)

            # ref) http://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
            if timeoutSec > 0:
                kill_proc = lambda proc: proc.kill()
                timer = Timer(timeoutSec, kill_proc, [p])

            # p.wait()

            try:
                stdOut, stdErr = p.communicate()
                print(stdOut)  # for printing slurm job id
                exitCode = p.returncode
            finally:
                if timeoutSec > 0:
                    timer.cancel()
                    exitCode = 143

            if runTime:
                end = time.time()
                elapsedSec = end - start
                if log:
                    log.info("*************************************")
                    if cmd.split(" ")[0].split("/")[-1]:
                        log.info(cmd.split(" ")[0].split("/")[-1])
                    log.info("Command took " + str(elapsedSec) + " sec.")

                    log.info("*************************************")

            if log and stdoutPrint:
                log.info(
                    "Return values: exitCode="
                    + str(exitCode)
                    + ", stdOut="
                    + str(stdOut)
                    + ", stdErr="
                    + str(stdErr)
                )

            if exitCode != 0:
                if log:
                    log.warn("- The exit code has non-zero value.")

    else:
        if log:
            log.error("- No command to run.")
            return None, None, -1

    return stdOut, stdErr, exitCode


# -------------------------------------------------------------------------------
def file_exist(path):
    """
    Check file existence
    """
    return os.path.exists(path)


# -------------------------------------------------------------------------------
def pad_string_path(myString, padLength=8, depth=None):
    """
    Pads a string with 0's and splits into groups of 2
    e.g. padStringPath(44) returns "00/00/00/44"

    @param myString: string to pad
    @return: padded string
    """
    myString = str(myString)
    padLength = int(padLength)

    if padLength > 8 or padLength <= 0:
        padLength = 8

    # left-pad with 0's
    myString = myString.zfill(padLength)

    # use re.findall function to split into pairs of strings
    stringList = re.findall("..", myString)

    # create ss/ss/ss/ss
    if not depth:
        padString = "/".join(stringList)
    else:
        padString = "/".join(stringList[:depth])

    padString = padString + "/"

    return padString
