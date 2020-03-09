#!/usr/bin/env python

"""
Command-line interface to Cromwell REST server.
This is used for testing Cromwell config; it is not part of JAWS.
"""

from datetime import datetime, timedelta
import sys
import os
import json
import requests
import click
import uuid
import re
import subprocess
import pathlib

if "TMPDIR" not in os.environ:
    sys.exit("TMPDIR env var not defined")
TMPDIR = os.environ["TMPDIR"]
if "CROMWELL_URL" not in os.environ:
    sys.exit("CROMWELL_URL env var not defined")
CROMWELL_URL = os.environ["CROMWELL_URL"]
if not CROMWELL_URL.startswith("http"):
    CROMWELL_URL = f"http://{CROMWELL_URL}"
WORKFLOWS_URL = f"{CROMWELL_URL}/api/workflows/v1"


@click.group()
def cromwell():
    """
    Cromwell REST server CLI
    """
    pass


def _get(url):
    try:
        r = requests.get(url)
    except Exception:
        sys.exit("Unable to GET from Cromwell at %s" % (url,))
    response = {}
    if r.status_code == 200:
        response["result"] = r.json()
    else:
        response["error"] = {"code": r.status_code, "message": "Run not found"}
    print(json.dumps(response, indent=4, sort_keys=True))
    return response


@cromwell.command()
@click.argument("run_id")
def status(run_id):
    """
    Retrieve the current status of a run.
    """
    url = "%s/%s/status" % (WORKFLOWS_URL, run_id)
    return _get(url)


@cromwell.command()
@click.argument("run_id")
def logs(run_id):
    """
    Retrieve the logs for a run.
    """
    url = "%s/%s/logs" % (WORKFLOWS_URL, run_id)
    return _get(url)


@cromwell.command()
@click.argument("run_id")
def outputs(run_id):
    """
    Retrieve the outputs for a run.
    """
    url = "%s/%s/outputs" % (WORKFLOWS_URL, run_id)
    return _get(url)


@cromwell.command()
@click.argument("run_id")
def metadata(run_id):
    """
    Retrieve the metadata of a run.
    """
    url = "%s/%s/metadata" % (WORKFLOWS_URL, run_id)
    return _get(url)


@cromwell.command()
@click.argument("run_id")
def cancel(run_id):
    """
    Cancel a run.
    """
    url = "%s/%s/abort" % (WORKFLOWS_URL, run_id)
    return _get(url)


@cromwell.command()
@click.argument("run_id")
def task_status(run_id):
    """
    Retrieve status of each task, with any error messages.
    """
    url = "%s/%s/metadata" % (WORKFLOWS_URL, run_id)
    r = requests.get(url)
    response = {}
    meta = None
    if r.status_code == 200:
        meta = r.json()
    else:
        response["error"] = {"code": r.status_code, "message": "Run not found"}
        print(json.dumps(response, indent=4, sort_keys=True))
        return response

    # load calls' status
    tasks = {}
    for task_name in meta["calls"]:
        start = ""
        try:
            start = meta["calls"][task_name][0]["start"]
        except Exception:
            start = "?"
        end = ""
        try:
            end = meta["calls"][task_name][0]["start"]
        except Exception:
            end = "?"
        tasks[start] = (
            task_name,
            meta["calls"][task_name][0]["executionStatus"],
            start,
            end,
        )

    # create sorted table
    result = []
    for start in reversed(sorted(tasks.keys())):
        result.append(tasks[start])

    response["result"] = result
    print(json.dumps(response, indent=4, sort_keys=True))
    return response


# simplifies the result JSON from that returned by cromwell
@cromwell.command()
@click.argument("run_id")
def get_labels(run_id):
    """
    Get labels for a run.
    """
    url = "%s/%s/labels" % (WORKFLOWS_URL, run_id)
    r = requests.get(url)
    response = {}
    if r.status_code == 200:
        result = r.json()
        labels = result["labels"]
        del labels["cromwell-workflow-id"]
        if "username" in labels:
            del labels["username"]
        response["result"] = labels
    else:
        response["error"] = {"code": 404, "message": "Run not found"}
    print(json.dumps(response, indent=4, sort_keys=True))
    return response


@cromwell.command()
@click.option("--username", default=os.environ["USER"], help="Username; default=USER")
def queue(username):
    """
    Get a user's unfinished runs.
    """
    status = ["Submitted", "Running", "Aborting"]
    labels = ["username:" + username]
    data = {"status": status, "label": labels}
    url = "%s/query" % (WORKFLOWS_URL,)
    r = requests.get(url, data)
    response = {}
    if r.status_code == 200:
        response["result"] = r.json()
    else:
        response["error"] = {"code": 404, "message": "Run not found"}
    print(json.dumps(response, indent=4, sort_keys=True))
    return response


@cromwell.command()
@click.option("--username", default=os.environ["USER"], help="Username; default=USER")
@click.option("--days", default=10, help="Delta days, default=10")
def history(username, days):
    """
    Get a user's run history.
    """
    d = datetime.today() - timedelta(int(days))
    start_date = d.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    labels = ["username:" + username]
    data = {"label": labels, "start": start_date}
    url = "%s/query" % (WORKFLOWS_URL,)
    r = requests.get(url, data)
    response = {}
    if r.status_code == 200:
        response["result"] = r.json()
    else:
        response["error"] = {"code": 404, "message": "Run not found"}
    print(json.dumps(response, indent=4, sort_keys=True))
    return response


@cromwell.command()
@click.argument("json_file")
@click.argument("wdl_file")
@click.option("--zip_file", default=None, help="subworkflows zip")
@click.option("--username", default=os.environ["USER"], help="Username; default=USER")
@click.option("--labels", default=None, help="JSON file of labels")
def submit(wdl_file, json_file, username, labels, zip_file):
    """
    Submit a run.
    """
    tmp_basename = os.path.join(TMPDIR, str(uuid.uuid4()))
    print("Using tmp base: %s" % (tmp_basename,))
    print("username : %s" % (username,))

    # labels
    tmp_labels_file = tmp_basename + ".labels.json"
    _prepare_labels_file(tmp_labels_file, username, labels)

    # SUBMIT TO CROMWELL
    files = {}
    files["workflowInputs"] = (
        "workflowInputs",
        open(json_file, "r"),
        "application/json",
    )
    files["workflowSource"] = (
        "workflowSource",
        open(wdl_file, "r"),
        "application/json",
    )
    files["labels"] = ("labels", open(tmp_labels_file, "r"), "application/json")
    if zip_file is not None:
        files["workflowDependencies"] = (
            "workflowDependencies",
            open(zip_file, "rb"),
            "application/zip",
        )

    print("Submitting run to %s" % (WORKFLOWS_URL,))
    r = requests.post(WORKFLOWS_URL, files=files)

    response = {}
    if r.status_code == 201:
        response["result"] = r.json()
    elif r.status_code == 500:
        response["error"] = {"code": r.status_code, "message": "Internal error"}
    elif r.status_code == 400:
        response["error"] = {
            "code": r.status_code,
            "message": "Invalid submission request",
        }
    else:
        response["error"] = {"code": r.status_code, "message": "Unknown error"}

    # CLEANUP
    os.remove(tmp_labels_file)

    print(json.dumps(response, indent=4, sort_keys=True))
    return response


def _prepare_labels_file(outfile, username, infile=None):
    """
    Given a username and optionally other labels (in a dict), create a JSON file.
    """
    assert outfile
    assert username
    labels = {}
    if infile is not None:
        with open(infile, "r") as f:
            data = json.load(f)
            for key in data:
                labels[key] = data[key]
    labels["username"] = username
    with open(outfile, "w") as f:
        json.dump(labels, f)


@cromwell.command()
@click.argument("src", nargs=1)
@click.argument("dest", nargs=1)
@click.option("--flatten", default=False)
def wfcopy(src, dest, flatten):
    """Reformat and copy cromwell run output"""
    src = os.path.abspath(src)
    dest = os.path.abspath(dest)
    shardname = None
    subwfname = None
    taskname = None
    cromwellFilesToSkip = [
        "stdout.background",
        "stderr.background",
        "script.background",
        "script.submit",
    ]

    logdir = os.path.join(dest, "log")
    if not os.path.exists(dest):
        os.makedirs(dest)
    if not os.path.exists(logdir):
        os.makedirs(logdir)

    rcfile = os.path.join(logdir, "workflow.rc")
    if os.path.exists(rcfile):
        os.remove(rcfile)
    with open(rcfile, "w") as fh:
        fh.write("#ExitCode\tTask\n")

    for rootDir, subdirs, files in os.walk(src):
        shardname = None
        parentDir = str(pathlib.Path(rootDir).parent)

        # look for call-* in directory name. If found, assign taskname.
        if os.path.basename(rootDir).startswith("call-") and rootDir == os.path.join(
            src, os.path.basename(rootDir)
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
            if not flatten and shardname:
                taskDir = os.path.join(dest, "%s/%s" % (taskname, shardname))
            else:
                taskDir = os.path.join(dest, taskname)

            if not os.path.exists(taskDir):
                os.makedirs(taskDir)

            for dname in subdirs:
                rsync(os.path.join(rootDir, dname), taskDir)

            for fname in files:
                fullname = os.path.join(rootDir, fname)
                outname = "%s-%s" % (taskname, shardname) if shardname else taskname

                if fname == "stdout":
                    rsync(fullname, os.path.join(logdir, "%s.stdout" % outname))
                if fname == "stderr":
                    rsync(fullname, os.path.join(logdir, "%s.stderr" % outname))
                if fname == "script":
                    rsync(fullname, os.path.join(logdir, "%s.script" % outname))
                if fname not in cromwellFilesToSkip:
                    rsync(fullname, taskDir)
                if fname == "rc":
                    with open(fullname, "r") as fh:
                        exitcode = fh.readlines()[0].strip()
                    with open(rcfile, "a") as fh:
                        fh.write("%s\t%s\n" % (exitcode, outname))


def rsync(src, dest):
    """ Copy source to destination using rsync

    :param src: Source path
    :type src: str
    :param dest" Destination path
    :type dest: str
    """
    cmd = "rsync -a %s %s" % (src, dest)
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    if process.returncode:
        sys.stderr.write(stderr.strip())
        sys.exit(process.returncode)


def getSubflowDirname(dirname):
    """
    Detect if input directory is a subworkflow (contains more than one 'call-' in the name).
    Ex: dirname=*/call-runMainBlast/blast.runblastplus/4a1eeaa3-79e7-42f0-9154-e5fb0636e0d8/
        call-runblastplus/shard-0  returns: runblastplus
    Ex: dirname=*/call-runMainBlast/blast.runblastplus/4a1eeaa3-79e7-42f0-9154-e5fb0636e0d8/
        call-runblastplus/execution  returns: runblastplus
    Ex: dirname=*/call-runMainBlast/execution  returns: None
    Ex: dirname=.../call-runMainBlast/shard-0  returns: None

    :param dirname: Path to subworkflow directory
    :type dirname: str
    :return: Subworkflow name if contains shards, None otherwise
    :rtype: str, None
    """
    subwfname = None
    if dirname.count("call-") > 1:
        rx = re.search(r"([^\/]+)\/shard-", dirname)
        subwfname = rx.group(1) if rx else os.path.basename(dirname)
        subwfname = subwfname.replace("call-", "", 1)
    return subwfname


if __name__ == "__main__":
    cromwell()
