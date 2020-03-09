"""
Utility to copy cromwell output to user specified output dir.
Input cromwell directory name must end with the job id
(e.g, cromwell-executions/fungal_min/4a14ec42-2d1c-4b65-9c4f-9f3a6307a094)
"""

import os
import re
import subprocess
import pathlib
import logging


logger = logging.getLogger(__package__)


def runCommand(cmd):
    """ Run command in shell. """


def rsync(src, dest):
    """ Copy source to destination using rsync.

    :param src: Source path which is Cromwell run output dir
    :type src: str
    :param dest: Destination path for reformatted output files
    :type dest: str
    :return:
    """
    cmd = "rsync -a %s %s" % (src, dest)
    _, stderr, exitcode = runCommand(cmd)
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    if process.returncode:
        logger.warn(f"Failed to rsync {src} to {dest}: " + stderr.strip())


def getSubflowDirname(dirname):
    """
    Detect if input directory is a subworkflow (contains more than one 'call-' in the name). If found, return
    subworkflow name:
    Ex: dirname=*/call-runMainBlast/blast.runblastplus/4a1eeaa3-79e7-42f0-9154-e5fb0636e0d8/
        call-runblastplus/shard-0 returns: runblastplus
    Ex: dirname=*/call-runMainBlast/blast.runblastplus/4a1eeaa3-79e7-42f0-9154-e5fb0636e0d8/
        call-runblastplus/execution returns: runblastplus
    Ex: dirname=*/call-runMainBlast/execution returns: None
    Ex: dirname=.../call-runMainBlast/shard-0 returns: None

    :param dirname: Path to Cromwell run output task folder
    :type dirname: str
    :return: Subworkflow name or None
    :rtype: str
    """
    subwfname = None
    if dirname.count("call-") > 1:
        rx = re.search(r"([^\/]+)\/shard-", dirname)
        subwfname = rx.group(1) if rx else os.path.basename(dirname)
        subwfname = subwfname.replace("call-", "", 1)
    return subwfname


def wfcopy(cromwellDir, outDir, flattenShardDir=False):
    """ Flatten cromwell folders and copy to new destination.

    :param cromwellDir: The base directory of the Cromwell run output.
    :type cromwellDir: str
    :param outDir: The destination directory of the reformatted output.
    :type outDir: str
    :param flattenShardDir: If True, shard output will be output to one dir, otherwise keep multiple subdirs.
    :type flattenShardDir: bool
    """
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
                rsync(os.path.join(rootDir, dname), taskDir)

            for fname in files:
                fullname = os.path.join(rootDir, fname)
                outname = "%s-%s" % (taskname, shardname) if shardname else taskname

                if fname == "stdout":
                    rsync(
                        fullname,
                        os.path.join(logdir, "%s.stdout" % outname),
                    )
                if fname == "stderr":
                    rsync(
                        fullname,
                        os.path.join(logdir, "%s.stderr" % outname),
                    )
                if fname == "script":
                    rsync(
                        fullname,
                        os.path.join(logdir, "%s.script" % outname),
                    )
                if fname not in cromwellFilesToSkip:
                    rsync(fullname, taskDir)
                if fname == "rc":
                    with open(fullname, "r") as fh:
                        exitcode = fh.readlines()[0].strip()
                    with open(rcfile, "a") as fh:
                        fh.write("%s\t%s\n" % (exitcode, outname))
