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


class OutputFolderExists(Exception):
    pass


def _rsync(src: str, dest: str) -> None:
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


def _get_subworkflow_name(dirname: str) -> str:
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


def _fix_perms(path: str) -> None:
    """Recursively chmod.

    :param path: Root dir
    :type dirname: str
    :return:
    """
    os.chmod(path, 0o0775)
    for dirpath, dirnames, filenames in os.walk(path):
        for dname in dirnames:
            os.chmod(os.path.join(dirpath, dname), 0o0775)
        for fname in filenames:
            os.chmod(os.path.join(dirpath, fname), 0o0664)


def get_files(srcdir: str, delimiter: str = "/", flatten_shard_dir: bool = False) -> tuple:
    """ Generator to traverse through a cromwell run directory. For each file in the cromwell task execution directory,
    return a tuple where the first element is the renamed task name and the second element is a file or directory
    within the cromwell task execution directory.

    :param srcdir: The base directory of the Cromwell run output.
    :type srcdir: str
    :param delimiter: the delimiter character used to concatenate the names of the cromwell task and subworkflow dirs.
    :type delimiter: str
    :param flatten_shard_dir: If True, shard output will be output to one dir, otherwise keep multiple subdirs.
    :type flatten_shard_dir: bool
    """

    for root_dir, subdirs, files in os.walk(srcdir):
        subtaskdir = None
        parent_dir = str(pathlib.Path(root_dir).parent)

        # if directory name starts with call-*, assign taskdir.
        if os.path.basename(root_dir).startswith("call-") and root_dir == os.path.join(
            srcdir, os.path.basename(root_dir)
        ):
            taskdir = re.sub(r"^call-", "", os.path.basename(root_dir))

        subwfname = _get_subworkflow_name(parent_dir)

        # if directory is a subworkflow, assign directory name with subworkflow name.
        if subwfname:
            # if subworkflow is a shard dir and flatten_shard_dir is False, append shard dir name to subworkflow name.
            if not flatten_shard_dir and os.path.basename(parent_dir).startswith("shard-"):
                subtaskdir = f"{subwfname}{delimiter}{os.path.basename(parent_dir)}"
            else:
                subtaskdir = subwfname

        # if not subworkflow but is a shard dir and flatten_shard_dir=False, assign directory name with shard name.
        elif os.path.basename(parent_dir).startswith("shard-") and not flatten_shard_dir:
            subtaskdir = os.path.basename(parent_dir)

        # if directory is 'execution', return task name and cromwell files within the dir.
        if os.path.basename(root_dir) == "execution":
            if subtaskdir:
                this_task = f"{taskdir}{delimiter}{subtaskdir}"
            else:
                this_task = taskdir
            for fname in files:
                yield this_task, os.path.join(os.path.abspath(root_dir), fname)
            for dname in subdirs:
                yield this_task, os.path.join(os.path.abspath(root_dir), dname)


def wfcopy(srcdir, dstdir, flatten_shard_dir=False):
    """ Given a cromwell run directory, copies all files and directories within the task execution dir to
    a renamed task directory in the new destination. The renamed directory is a flattened representation of the
    cromwell tasks including subworkflows and shard directories.

    :param srcdir: The base directory of the Cromwell run output.
    :type srcdir: str
    :param dstdir: The destination directory of the reformatted output.
    :type dstdir: str
    :param flatten_shard_dir: If True, shard output will be output to one dir, otherwise keep multiple subdirs.
    :type flatten_shard_dir: bool
    """
    if os.path.exists(dstdir):
        raise OutputFolderExists(dstdir)

    # these are cromwell files to store in output path's log dir.
    target_log_files = [
        "script",
        "script.submit",
        "stdout",
        "stdout.submit",
        "stderr",
        "stderr.submit",
    ]
    logdir = os.path.join(dstdir, "log")
    logname_delimiter = "."
    taskname_delimiter = "/"
    rcfile = os.path.join(logdir, "workflow.rc")

    if not os.path.exists(dstdir):
        os.makedirs(dstdir)
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    if os.path.exists(rcfile):
        os.remove(rcfile)
    with open(rcfile, "w") as fh:
        fh.write("#ExitCode\tTask\n")

    # get all files/dirs from the cromwell execution directory along with the wfcopy formatted task name.
    for taskname, filename in get_files(
        srcdir, delimiter=taskname_delimiter, flatten_shard_dir=flatten_shard_dir
    ):
        dst_taskdir = os.path.join(dstdir, taskname)
        basename = os.path.basename(filename)

        if not os.path.exists(dst_taskdir):
            os.makedirs(dst_taskdir)

        # if filename is a cromwell file for logging, rename and copy file to log directory.
        if basename in target_log_files and os.path.isfile(filename):
            dst_taskdir = os.path.join(
                logdir,
                f"{taskname.replace(taskname_delimiter, logname_delimiter)}.{basename}",
            )

        # if filename is 'rc', extract return code and write to log directory's rc file.
        if basename == "rc":
            with open(filename, "r") as fh:
                exitcode = fh.readlines()[0].strip()
            with open(rcfile, "a") as fh:
                fh.write(
                    "%s\t%s\n"
                    % (
                        exitcode,
                        taskname.replace(taskname_delimiter, logname_delimiter),
                    )
                )
        else:
            _rsync(filename, dst_taskdir)

    # change output directory permissions.
    _fix_perms(dstdir)
