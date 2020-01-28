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
Created on Apr 3, 2013

Convenience methods to make it easier to run external programs
and other os-related tools

@author: Seung-Jin Sul

"""

from subprocess import Popen, call, PIPE
import os
import glob
import sys
import shlex
import unittest
import errno
import re

defineCalledProcessError = False
try:
    from subprocess import CalledProcessError
except ImportError:
    defineCalledProcessError = True

if defineCalledProcessError:
    class CalledProcessError(OSError):
        def __init__(self, returncode, cmd, *l, **kw):
            OSError.__init__(self, *l, **kw)
            self.cmd = cmd
            self.returncode = returncode


#-------------------------------------------------------------------------------
def run(*popenargs, **kwargs):
    """
    Run user command using subprocess.call
    :param popenargs: command and options to run
    :param kwargs: additional parameters
    :return:
    """
    kw = {}
    kw.update(kwargs)
    dryRun = kw.pop('dryRun', False)

    if dryRun:
        print(popenargs)
    else:
        ## convert something like run("ls -l") into run("ls -l", shell=True)
        if isinstance(popenargs[0], str) and len(shlex.split(popenargs[0])) > 1:
            kw.setdefault("shell", True)

        ## > /dev/null 2>&1
        if kw.pop("supressAllOutput", False):
            stdnull = open(os.devnull, "w") ## incompat with close_fds on Windows
            kw.setdefault("stdout", stdnull)
            kw.setdefault("stderr", stdnull)
        else:
            stdnull = None

        returncode = call(*popenargs, **kw)
        if stdnull:
            stdnull.close()
        if returncode != 0:
            raise CalledProcessError(returncode=returncode, cmd=str(popenargs))


#-------------------------------------------------------------------------------
def backTicks(*popenargs, **kwargs):
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
    dryRun = kw.pop('dryRun', False)
    dryRet = kw.pop('dryRet', None)

    if dryRun:
        print(popenargs)
        return dryRet
    else:
        kw['stdout'] = PIPE
        p = Popen(*popenargs, **kw)
        retOut = p.communicate()[0]
        if p.returncode != 0:
            raise CalledProcessError(returncode=p.returncode, cmd=str(popenargs))
        return retOut


#-------------------------------------------------------------------------------
def makeDir(path, dryRun=False):
    """
    Create one dir with pathname path or do nothing if it already exists. Same as Linux 'mkdir -p'.
    :param path: path to create
    :param dryRun: dryrun directive
    :return:
    """
    if not dryRun:
        if not os.path.exists(path):
            os.makedirs(path)
    else:
        print("makeDir %s" % (path, ))


def make_dir_p(path):
    """
    The method make_dir_p() is recursive directory creation function.
    Like mkdir(), but makes all intermediate-level directories needed to contain the leaf directory.
    :param path: path to create
    :return:
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


#-------------------------------------------------------------------------------
def makeDirs(paths, dryRun=False):
    """
    Create muiltiple dirs with the same semantics as makeDir
    :param paths: path to delete
    :param dryRun: dryrun directive
    :return:
    """
    for path in paths:
        makeDir(path=path, dryRun=dryRun)


#-------------------------------------------------------------------------------
def makeFilePath(fileName):
    """
    Assume that the argument is a file name and make all directories that are
    part of it
    :param fileName: create dir to the file
    :return:
    """
    dirName = os.path.dirname(fileName)
    if dirName not in ("", "."):
        makeDir(dirName)


#-------------------------------------------------------------------------------
def rmDir(path, dryRun=False):
    """
    Remove dir
    :param path: path to delete
    :param dryRun: dryrun directive
    :return:
    """
    ## To do: perhaps use shutil.rmtree instead?
    run(["rm", "-rf", path], dryRun=dryRun)

## make alias
rmrf = rmDir


#-------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------
def rmfMany(paths, dryRun=False):
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


#-------------------------------------------------------------------------------
def remakeDir(path, dryRun=False):
    """
    Create an empty dir with a given path.
    If path already exists,  it will be removed first.
    :param path: path to delete and create
    :param dryRun:
    :return:
    """
    rmrf(path, dryRun=dryRun)
    makeDir(path, dryRun=dryRun)


#-------------------------------------------------------------------------------
def chmod(path, mode, opts='', dryRun=False):
    """
    Change mode.
    :param path: path to chmod
    :param mode: the form `[ugoa]*([-+=]([rwxXst]*|[ugo]))+'.
    :param opts: additional chmod options
    :param dryRun: dryrun directive
    :return:
    """
    if isinstance(path, basestring):
        path = [path]
    else:
        path = list(path)
    run(["chmod"]+opts.split()+[mode]+path, dryRun=dryRun)


#-------------------------------------------------------------------------------
def run_sh_command(cmd, live=False, log=None, runTime=False, stdoutPrint=True, timeoutSec=0):
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
    jobid = -1

    if cmd:
        if not live:
            stdOut = "Not live: cmd = '%s'" % (cmd)
            exitCode = 0
        else:
            if log and stdoutPrint:
                log.info("cmd: %s" % (cmd))

            ##---------
            ## OLD
            if runTime:
                start = time.time()

            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)

            ## ref) http://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
            if timeoutSec > 0:
                kill_proc = lambda proc: proc.kill()
                timer = Timer(timeoutSec, kill_proc, [p])

            #p.wait()

            try:
                stdOut, stdErr = p.communicate()
                print(stdOut)  # for printing slurm job id
                exitCode = p.returncode
            finally:
                if timeoutSec > 0:
                    timer.cancel()
                    exitCode = 143
                else:
                    pass

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
                log.info("Return values: exitCode=" + str(exitCode) + ", stdOut=" + str(stdOut) + ", stdErr=" + str(stdErr))

            if exitCode != 0:
                if log:
                    log.warn("- The exit code has non-zero value.")

    else:
        if log:
            log.error("- No command to run.")
            return None, None, -1

    return stdOut, stdErr, exitCode


#-------------------------------------------------------------------------------
def fExist(path, dryRun=False):
    """
    Check file existence
    """
    return os.path.exists(path)


#-------------------------------------------------------------------------------
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


# # -------------------------------------------------------------------------------
# def jamo_fetch_cmd(metadata_id):
#     fetch_cmd = ""
#     file_id = ""
#     file_name = ""
#
#     if metadata_id:
#         sdm_url = "https://sdm2.jgi-psf.org"
#         sdm_api = "/api/metadata/file/[metadata_id]"
#         sdm_api = sdm_api.replace("[metadata_id]", metadata_id)
#
#         curl = Curl(sdm_url)
#         response = None
#
#         try:
#             response = curl.get(sdm_api)
#         except CurlHttpException as e:
#             print e.response
#             # transLog.error("- Failed SDM API Call: %s/%s, %s  %s", sdmURL, sdmRawFileAPI, e, RQCGlobals.msgFail)
#         except Exception as e:
#             print e.args
#
#         if response:
#             if type(response) == str:
#                 response = json.loads(response)
#             # sys.exit(44)
#             if 'file_id' in response:
#                 # print "****** %s" % type(response['file_id']) # weird response in denovo
#                 fetch_cmd = "jamo fetch custom file_id=%s" % (str(response['file_id']))
#                 file_id = response['file_id']
#
#             if 'current_location' in response:
#                 file_name = response['current_location']
#
#     fetch_cmd += " # restore to: %s" % file_name
#
#
#     return fetch_cmd


#-------------------------------------------------------------------------------
class TestOsUtility(unittest.TestCase):
    def testRun(self):
        try:
            run(["rm", "-rf", "./unittest"], dryRun=False)
        except CalledProcessError as msg:
            self.assertNotEqual(msg.returncode, 0)
        try:
            makeDir("./unittest", dryRun=False)
        except CalledProcessError as msg:
            self.assertEqual(msg.returncode, 0)
        try:
            rmDir("./unittest", dryRun=False)
        except CalledProcessError as msg:
            self.assertEqual(msg.returncode, 0)

    def testBackTicks(self):
        cmd = "free"
        try:
            freeOut = backTicks(cmd, shell=True)
        except CalledProcessError as msg:
            sys.stderr.write("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
            sys.exit(1)

        ret = -1
        ret = float(freeOut.split('\n')[1].split()[2]) / \
              float(freeOut.split('\n')[1].split()[1]) * 100.0
        assert ret > -1


#
# Use this for testing purpose
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    # unittest.main()