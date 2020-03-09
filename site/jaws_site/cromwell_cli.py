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

TMPDIR = os.environ["TMPDIR"]
CROMWELL_PORT = os.environ["CROMWELL_PORT"]
WORKFLOWS_URL = f"http://localhost:{CROMWELL_PORT}/api/workflows/v1"


@click.group()
def cromwell():
    """Cromwell REST server CLI"""
    pass


def _get(url):
    """REST "GET" function helper

    :param url: URL to GET
    :type url: str
    :return: response from server, decoded from JSON string
    :rtype: dict
    """
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
    """Retrieve the current status of a run.

    :param run_id: Cromwell run ID
    :type run_id: str
    :return: Response JSON includes status key.
    :rtype: dict
    """
    url = "%s/%s/status" % (WORKFLOWS_URL, run_id)
    return _get(url)


@cromwell.command()
@click.argument("run_id")
def logs(run_id):
    """Retrieve the logs for a run.

    :param run_id: Cromwell run ID
    :type run_id: str
    :return: Cromwell run log
    :rtype: str
    """
    url = "%s/%s/logs" % (WORKFLOWS_URL, run_id)
    return _get(url)


@cromwell.command()
@click.argument("run_id")
def outputs(run_id):
    """Retrieve the outputs for a run.

    :param run_id: Cromwell run ID
    :type run_id: str
    :return: Dictionary indicating paths to output files created by the run
    :rtype: dict
    """
    url = "%s/%s/outputs" % (WORKFLOWS_URL, run_id)
    return _get(url)


@cromwell.command()
@click.argument("run_id")
def metadata(run_id):
    """Retrieve the metadata of a run.

    :param run_id: Cromwell run ID
    :type run_id: str
    :return: Cromwell metadata for a run, if found.
    :rtype: dict
    """
    url = "%s/%s/metadata" % (WORKFLOWS_URL, run_id)
    return _get(url)


@cromwell.command()
@click.argument("run_id")
def cancel(run_id):
    """Cancel a run.

    :param run_id: Cromwell run ID
    :type run_id: str
    :return: Cromwell's response, indicating success or failure.
    :rtype: dict
    """
    url = "%s/%s/abort" % (WORKFLOWS_URL, run_id)
    return _get(url)


@cromwell.command()
@click.argument("run_id")
def task_status(run_id):
    """Retrieve status of each task, with any error messages.

    :param run_id: Cromwell run ID
    :type run_id: str
    :return: Dictionary of task and status, with any error messages.
    :rtype: dict
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
    """Get labels for a run.

    :param run_id: Cromwell run ID
    :type run_id: str
    :return: Dictionary of labels associated with a Cromwell run.
    :rtype: dict
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
    """Get a user's unfinished runs.

    :param username: Value of label username to search.
    :type run_id: str
    :return: List of unfinished runs owned by specified user.
    :rtype: dict
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
    """Get a user's run history.

    :param username: Value of label username to search.
    :type run_id: str
    :param days: The time window in which to search, from now.
    :type days: int
    :return: a list of runs within a time window, belonging to specified user
    :rtype: dict
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
    """Submit a run.

    :param wdl_file: Path to WDL file for the run
    :type wdl_file: str
    :param json_file: Path to inputs-JSON file for the run
    :type json_file: str
    :param username: Value of label username to associate with run.
    :type run_id: str
    :param labels: Path to JSON file containing dict of label kv-pairs to associate with run.
    :type labels: str
    :param zip_file: A zip file of all subworkflows required by the WDL.
    :type zip_file: str
    :return: Response from Cromwell
    :rtype: dict
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
    """Given a username and optionally other labels (in a dict), create a JSON file.

    :param outfile: Path to JSON-formatted output file of labels to associate with the run.
    :type outfile: str
    :param username: Name of user to associate with the run, using username label.
    :type username: str
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


if __name__ == "__main__":
    cromwell()
