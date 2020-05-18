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

from jaws_site import wfcopy as wfc

if "USER" not in os.environ:
    sys.exit("USER env var not defined")
if "TMPDIR" not in os.environ:
    sys.exit("TMPDIR env var not defined")
TMPDIR = os.environ["TMPDIR"]
if "CROMWELL_URL" not in os.environ:
    sys.exit("CROMWELL_URL env var not defined")
CROMWELL_URL = os.environ["CROMWELL_URL"]
if not CROMWELL_URL.startswith("http"):
    CROMWELL_URL = f"http://{CROMWELL_URL}"
WORKFLOWS_URL = f"{CROMWELL_URL}/api/workflows/v1"


def _prepare_labels_file(outfile, username, infile=None):
    """Given a username and optional dict, create a JSON file."""
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


def _post(url, data={}, files={}, tmpfiles=[]):
    """POST to URL, print result, exit on error."""
    try:
        r = requests.post(url, data=data, files=files)
    except requests.exceptions.RequestException as e:
        raise e(f"Unable to POST to {url}: {e}")
    for tmpfile in tmpfiles:
        os.remove(tmpfile)
    if r.status_code == 201:
        print(json.dumps(r.json(), indent=4, sort_keys=True))
        return r.json()
    else:
        sys.exit(r.text)


def _get(url, data={}):
    """GET from URL, print result, exit on error."""
    try:
        r = requests.get(url, data)
    except requests.exceptions.RequestException as e:
        raise e(f"Unable to GET from {url}: {e}")
    if r.status_code == 200:
        print(json.dumps(r.json(), indent=4, sort_keys=True))
        return r.json()
    else:
        sys.exit(r.text)


@click.group()
def cromwell():
    """Cromwell REST server CLI"""
    pass


@cromwell.command()
@click.argument("run_id")
def status(run_id):
    """Retrieve the current status of a run."""
    _get(f"{WORKFLOWS_URL}/{run_id}/status")


@cromwell.command()
@click.argument("run_id")
def logs(run_id):
    """Retrieve the logs for a run."""
    _get(f"{WORKFLOWS_URL}/{run_id}/logs")


@cromwell.command()
@click.argument("run_id")
def outputs(run_id):
    """Retrieve the outputs for a run."""
    _get(f"{WORKFLOWS_URL}/{run_id}/outputs")


@cromwell.command()
@click.argument("run_id")
def metadata(run_id):
    """Retrieve a run's metadata."""
    _get(f"{WORKFLOWS_URL}/{run_id}/metadata")


@cromwell.command()
@click.argument("run_id")
def abort(run_id):
    """Abort a run."""
    _post(f"{WORKFLOWS_URL}/{run_id}/abort")


@cromwell.command()
@click.argument("run_id")
def task_status(run_id):
    """Retrieve status of each task."""
    url = "%s/%s/metadata" % (WORKFLOWS_URL, run_id)
    r = requests.get(url)
    meta = None
    if r.status_code == 200:
        meta = r.json()
    else:
        sys.exit(r.text)
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
    result = []
    for start in reversed(sorted(tasks.keys())):
        result.append(tasks[start])
    print(json.dumps(result, indent=4, sort_keys=True))


@cromwell.command()
@click.argument("run_id")
def get_labels(run_id):
    """Get labels for a run."""
    _get(f"{WORKFLOWS_URL}/{run_id}/labels")


@cromwell.command()
@click.option("--username", default=os.environ["USER"], help="Username; default=USER")
def queue(username):
    """Get a user's unfinished runs."""
    data = {
        "status": ["Submitted", "Running", "Aborting"],
        "label": ["username:" + username],
    }
    _get(f"{WORKFLOWS_URL}/query", data=data)


@cromwell.command()
@click.option("--username", default=os.environ["USER"], help="Username; default=USER")
@click.option("--days", default=10, help="Delta days, default=10")
def history(username, days):
    """Get a user's run history."""
    d = datetime.today() - timedelta(int(days))
    data = {
        "label": ["username:" + username],
        "start": d.strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
    }
    _get(f"{WORKFLOWS_URL}/query", data=data)


@cromwell.command()
@click.argument("json_file")
@click.argument("wdl_file")
@click.option("--zip_file", default=None, help="subworkflows zip")
@click.option("--username", default=os.environ["USER"], help="Username; default=USER")
@click.option("--labels", default=None, help="JSON file of labels")
def submit(wdl_file, json_file, username, labels, zip_file):
    """Submit a run."""
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
    _post(WORKFLOWS_URL, files=files, tmpfiles=[tmp_labels_file])


@cromwell.command()
@click.argument("src", nargs=1)
@click.argument("dest", nargs=1)
@click.option("-f", "--flatten", is_flag=True)
def wfcopy(src, dest, flatten):
    """Reformat and copy cromwell run output"""
    src = os.path.abspath(src)
    dest = os.path.abspath(dest)
    wfc.wfcopy(src, dest, flatten_shard_dir=flatten)
