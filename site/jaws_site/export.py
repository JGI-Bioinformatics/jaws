#!/usr/bin/env python

"""
This script periodically checks the workflow manager (Cromwell) for newly completed runs (i.e. not yet exported).  The output files are rearranged into concise format, packaged into a tarball, and sent to the appropriate Globus endpoint.  If the run was successful, the tmpfiles are also purged.  The state of each run is tracked in a MySQL database.
"""

import sys
import os
import db

def get_incomplete_runs()
    """
    Retrieve a list of runs to watch out for, from the MySQL database.
    """

def get_run_status():
    """
    Ask Cromwell for the status of a run.
    """

def query_workflow_manager():
    """
    Ask Cromwell for the status of each watched run.
    """
    for run_id in runs:
        status = get_run_status(run_id)
        if status == 'COMPLETED':
            foo
        elif status == 'CANCELED':
            bar
        elif status == 'FAILED':
            fubar
        elif status == 'RUNNING':
            1 # do nothing
        else:
            sys.stderr.write("Unrecognized job status: %s" % (status,))

def get_transfers_in_progress():
    """
    Retrieve a list of unfinished Globus task ids from the MySQL database.
    """
    stmt = 'SELECT globus_task_id FROM transfers WHERE state = "DOWNLOADING"'

def get_transfer_status():
    """
    Ask Globus for the state of one file transfer.
    """

def query_globus():
    """
    Ask Globus for the state of each file transfer.
    """
    for globus_task_id in uploads:
        status = get_transfer_status(globus_task_id)

def package_outputs(run_dir, tar_file)
    """
    Run outputs are reformatted and packaged into a tarball.
    """

def xfer_tarball(tar_file, dest_endpoint):
    """
    Transfer the tarball to the indicated Globus endpoint.  Block until done.
    """
    # TODO
    return update_run_state(run_id, "READY")

def update_run_state(run_id, new_state):
    """
    Change the run's state.
    """
    stmt = "UPDATE run SET state = %s WHERE run_id = %s"
    try:
        # TODO
    except:
        # TODO

@util.command()
@click.argument("cromwelldir")
@click.argument("outdir")
def wfcopy(cromwelldir, outdir, verbose=False):
    """
    Copy cromwell output to user specified dir.
    """
    logdir = os.path.join(outdir, "log")

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(logdir):
        os.makedirs(logdir)

    cromwellFilesToSkip = ['stdout.background', 'stderr.background', 'script.background', 'script.submit']
    taskname = None
    rcfile = os.path.join(logdir, "workflow.rc")

    if os.path.exists(rcfile):
        os.remove(rcfile)
    with open(rcfile, 'w') as fh:
        fh.write("#ExitCode\tTask\n")

    for rootdir, subdirs, files in os.walk(cromwelldir):
        if os.path.basename(rootdir).startswith('call-') and rootdir == os.path.join(cromwelldir, os.path.basename(rootdir)):
            taskname = re.sub(r'^call-', '', os.path.basename(rootdir))

        if rootdir.endswith('execution'):
            parentdir = str(pathlib.Path(rootdir).parent)
            shardname = None

            if re.search(r'shard-', parentdir):
                shardname = os.path.basename(parentdir)

            taskdir = os.path.join(outdir, taskname)
            if not os.path.exists(taskdir):
                os.makedirs(taskdir)

            for dname in subdirs:
                rsync(os.path.join(rootdir, dname), taskdir, verbose=verbose)

            for fname in files:
                fullname = os.path.join(rootdir, fname)
                outname = "%s-%s"%(taskname, shardname) if shardname else taskname

                if fname == 'stdout':
                    rsync(fullname, os.path.join(logdir, "%s.stdout"%outname), verbose=verbose)
                elif fname == 'stderr':
                    rsync(fullname, os.path.join(logdir, "%s.stderr"%outname), verbose=verbose)
                elif fname == 'script':
                    rsync(fullname, os.path.join(logdir, "%s.script"%outname), verbose=verbose)
                elif fname == 'rc':
                    with open(fullname, 'r') as fh:
                        exitcode = fh.readlines()[0].strip()
                    with open(rcfile, 'a') as fh:
                        fh.write("%s\t%s\n"%(exitcode, outname))
                elif fname not in cromwellFilesToSkip:
                    rsync(fullname, taskdir, verbose=verbose)

def runCommand(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.strip(), stderr.strip(), process.returncode

def rsync(src, dest, verbose=False):
    cmd = "rsync -a %s %s"%(src, dest)
    if verbose:
        print(cmd)
    stdout, stderr, exitcode = runCommand(cmd)
    if exitcode:
        sys.stderr.write(stderr)
        sys.exit(exitcode)
