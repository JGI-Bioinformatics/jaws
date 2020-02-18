#!/usr/bin/env python
"""
Utility to copy cromwell output to user specified output dir.
Input cromwell directory name must end with the job id (e.g,
cromwell-executions/fungal_min/4a14ec42-2d1c-4b65-9c4f-9f3a6307a094)

20190315 - initial version
20190805 - when calling a sub-workflow, use the root's task name instead of the sub workflow's task name.
20200121 - update code to preserve shard directory.
"""

import os
import sys
import re
import subprocess
import pathlib
import argparse

__version__ = "1.0.0"


def getArgs():
    """
    Parse command line arguments.
    """
    progDesc = """
    Utility to copy cromwell outputs to user specified output dir.
    """
    parser = argparse.ArgumentParser(
        description=progDesc, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-fs",
        "--flattenShardDir",
        dest="flattenShardDir",
        action="store_true",
        help="copy contents of shard dirs into one dir.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="print messages to stdout.",
    )
    parser.add_argument(
        "-V", "--version", action="version", version=__version__, help="show version"
    )
    parser.add_argument("srcpath", help="cromwell jobid path")
    parser.add_argument("dstpath", help="destination path")

    # If no input arguments specified, then display help menu.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


def runCommand(cmd):
    """ Run command in shell. """
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    return stdout.strip(), stderr.strip(), process.returncode


def rsync(src, dest, verbose=False):
    """ Copy source to destination using rsync. """
    cmd = "rsync -a %s %s" % (src, dest)
    if verbose:
        print(cmd)
    _, stderr, exitcode = runCommand(cmd)
    if exitcode:
        sys.stderr.write(stderr)
        sys.exit(exitcode)


def getSubflowDirname(dirname):
    """
    Detect if input directory is a subworkflow (contains more than one 'call-' in the name). If found, return
    subworkflow name:
    Ex: dirname=*/call-runMainBlast/blast.runblastplus/4a1eeaa3-79e7-42f0-9154-e5fb0636e0d8/call-runblastplus/shard-0
        returns: runblastplus

        dirname=*/call-runMainBlast/blast.runblastplus/4a1eeaa3-79e7-42f0-9154-e5fb0636e0d8/call-runblastplus/execution
        returns: runblastplus

        dirname=*/call-runMainBlast/execution
        returns: None

        dirname=.../call-runMainBlast/shard-0
        returns: None
    """
    subwfname = None
    if dirname.count("call-") > 1:
        rx = re.search(r"([^\/]+)\/shard-", dirname)
        subwfname = rx.group(1) if rx else os.path.basename(dirname)
        subwfname = subwfname.replace("call-", "", 1)
    return subwfname


def wfcopy(cromwellDir, outDir, flattenShardDir=False, verbose=False):
    """ Flatten cromwell folders and copy to new destination. """
    shardname = None
    subwfname = None
    taskname = None
    cromwellFilesToSkip = [
        "stdout.background",
        "stderr.background",
        "script.background",
        "script.submit",
    ]

    logdir = os.path.join(outDir, "log")
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    if not os.path.exists(logdir):
        os.makedirs(logdir)

    rcfile = os.path.join(logdir, "workflow.rc")
    if os.path.exists(rcfile):
        os.remove(rcfile)
    with open(rcfile, "w") as fh:
        fh.write("#ExitCode\tTask\n")

    for rootDir, subdirs, files in os.walk(cromwellDir):
        shardname = None
        parentDir = str(pathlib.Path(rootDir).parent)

        # look for call-* in directory name. If found, assign taskname.
        if os.path.basename(rootDir).startswith("call-") and rootDir == os.path.join(
            cromwellDir, os.path.basename(rootDir)
        ):
            taskname = re.sub(r"^call-", "", os.path.basename(rootDir))

        # look for shard directory. If found, assign shardname.
        if os.path.basename(parentDir).startswith("shard-"):
            shardname = os.path.basename(parentDir)
            subwfname = getSubflowDirname(parentDir)
            if subwfname:
                shardname = "%s-%s" % (subwfname, shardname)

        # look for execution directory. If found, copy files to destination.
        if rootDir.endswith("execution"):
            if not flattenShardDir and shardname:
                taskDir = os.path.join(outDir, "%s/%s" % (taskname, shardname))
            else:
                taskDir = os.path.join(outDir, taskname)

            if not os.path.exists(taskDir):
                os.makedirs(taskDir)

            for dname in subdirs:
                rsync(os.path.join(rootDir, dname), taskDir, verbose=verbose)

            for fname in files:
                fullname = os.path.join(rootDir, fname)
                outname = "%s-%s" % (taskname, shardname) if shardname else taskname

                if fname == "stdout":
                    rsync(
                        fullname,
                        os.path.join(logdir, "%s.stdout" % outname),
                        verbose=verbose,
                    )
                if fname == "stderr":
                    rsync(
                        fullname,
                        os.path.join(logdir, "%s.stderr" % outname),
                        verbose=verbose,
                    )
                if fname == "script":
                    rsync(
                        fullname,
                        os.path.join(logdir, "%s.script" % outname),
                        verbose=verbose,
                    )
                if fname not in cromwellFilesToSkip:
                    rsync(fullname, taskDir, verbose=verbose)
                if fname == "rc":
                    with open(fullname, "r") as fh:
                        exitcode = fh.readlines()[0].strip()
                    with open(rcfile, "a") as fh:
                        fh.write("%s\t%s\n" % (exitcode, outname))


if __name__ == "__main__":
    argParser = getArgs()
    cromwellDir = os.path.abspath(argParser.srcpath)
    outDir = os.path.abspath(argParser.dstpath)
    wfcopy(
        cromwellDir,
        outDir,
        flattenShardDir=argParser.flattenShardDir,
        verbose=argParser.verbose,
    )
    sys.exit()
