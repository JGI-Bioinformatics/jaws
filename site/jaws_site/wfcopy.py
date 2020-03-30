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


def rsync(src, dest):
    """ Copy source to destination using rsync.

    :param src: Source path which is Cromwell run output dir
    :type src: str
    :param dest: Destination path for reformatted output files
    :type dest: str
    :return:
    """
    cmd = "rsync -a %s %s" % (src, dest)
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


def wfcopy(in_dir, out_dir, flattenShardDir=False):
    """ Flatten cromwell folders and copy to new destination.

    :param in_dir: The base directory of the Cromwell run output.
    :type in_dir: str
    :param out_dir: The destination directory of the reformatted output.
    :type out_dir: str
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

    log_dir = os.path.join(out_dir, "log")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    rcfile = os.path.join(log_dir, "workflow.rc")
    if os.path.exists(rcfile):
        os.remove(rcfile)
    with open(rcfile, "w") as fh:
        fh.write("#ExitCode\tTask\n")

    for root_dir, subdirs, files in os.walk(in_dir):
        shardname = None
        parent_dir = str(pathlib.Path(root_dir).parent)

        # look for call-* in directory name. If found, assign taskname.
        if os.path.basename(root_dir).startswith("call-") and root_dir == os.path.join(
            in_dir, os.path.basename(root_dir)
        ):
            taskname = re.sub(r"^call-", "", os.path.basename(root_dir))

        # look for shard directory. If found, assign shardname.
        if os.path.basename(parent_dir).startswith("shard-"):
            shardname = os.path.basename(parent_dir)
            subwfname = getSubflowDirname(parent_dir)
            if subwfname:
                shardname = "%s-%s" % (subwfname, shardname)

        # look for execution directory. If found, copy files to destination.
        if root_dir.endswith("execution"):
            if not flattenShardDir and shardname:
                task_dir = os.path.join(out_dir, "%s/%s" % (taskname, shardname))
            else:
                task_dir = os.path.join(out_dir, taskname)

            if not os.path.exists(task_dir):
                os.makedirs(task_dir)

            for dname in subdirs:
                rsync(os.path.join(root_dir, dname), task_dir)

            for fname in files:
                fullname = os.path.join(root_dir, fname)
                outname = "%s-%s" % (taskname, shardname) if shardname else taskname

                if fname == "stdout":
                    rsync(
                        fullname, os.path.join(log_dir, "%s.stdout" % outname),
                    )
                if fname == "stderr":
                    rsync(
                        fullname, os.path.join(log_dir, "%s.stderr" % outname),
                    )
                if fname == "script":
                    rsync(
                        fullname, os.path.join(log_dir, "%s.script" % outname),
                    )
                if fname not in cromwellFilesToSkip:
                    rsync(fullname, task_dir)
                if fname == "rc":
                    with open(fullname, "r") as fh:
                        exitcode = fh.readlines()[0].strip()
                    with open(rcfile, "a") as fh:
                        fh.write("%s\t%s\n" % (exitcode, outname))
