#!/usr/bin/env python

"""
Extended Python API and CLI for Cromwell.  Primarily adds user support and convenience functions, such as zipping subworkflows.  Note that this package does not provide any authentication.  This is used only for testing.
"""

from datetime import datetime, timedelta
import sys
import os
import json
import requests
import click
import uuid
import zipfile

from jaws_site import config

verbose=True

config = config.Config(env="JAWS_SITE_CONFIG")

@click.group()
def jaws():
    """
    JAWS Site CLI
    """
    pass

def _get(url):
    r = requests.get(url)
    response = {}
    if r.status_code == 200:
        response["result"] = r.json()
    else:
        response["error"] = {
            "code" : r.status_code,
            "message" : "Run not found"
        }
    if verbose:
        print(json.dumps(response, indent=4, sort_keys=True))
    return response

@jaws.command()
@click.argument('run_id')
def status(run_id):
    """
    Retrieve the current status of a run.
    """
    url = config["CROMWELL"]['url'] + '/' + run_id + "/status"
    return _get(url)

@jaws.command()
@click.argument('run_id')
def logs(run_id):
    """
    Retrieve the logs for a run.
    """
    url = config["CROMWELL"]['url'] + '/' + run_id + "/logs"
    return _get(url)

@jaws.command()
@click.argument('run_id')
def outputs(run_id):
    """
    Retrieve the outputs for a run.
    """
    url = config["CROMWELL"]['url'] + '/' + run_id + "/outputs"
    return _get(url)

@jaws.command()
@click.argument('run_id')
def metadata(run_id):
    """
    Retrieve the metadata of a run.
    """
    url = config["CROMWELL"]['url'] + '/' + run_id + "/metadata"
    return _get(url)

@jaws.command()
@click.argument('run_id')
def cancel(run_id):
    """
    Cancel a run.
    """
    url = config["CROMWELL"]['url'] + '/' + run_id + "/abort"
    return _get(url)

@jaws.command()
@click.argument('run_id')
def get_labels(run_id):
    """
    Retrieve labels for a run.
    """
    url = config["CROMWELL"]['url'] + '/' + run_id + '/labels'
    return _get(url)

@jaws.command()
@click.argument('run_id')
def task_status(run_id):
    """
    Retrieve status of each task, with any error messages.
    """
    url = config["CROMWELL"]['url'] + '/' + run_id + "/metadata"
    r = requests.get(url)
    response = {}
    meta = None
    if r.status_code == 200:
        meta = r.json()
    else:
        response["error"] = {
            "code" : r.status_code,
            "message" : "Run not found"
        }
        if verbose:
            print(json.dumps(response, indent=4, sort_keys=True))
        return response

    # load calls' status
    tasks = {}
    for task_name in meta['calls']:
        start = ''
        try:
            start = meta['calls'][task_name][0]['start']
        except:
            start = '?'
        end = ''
        try:
            end = meta['calls'][task_name][0]['start']
        except:
            end = '?'
        tasks[start] = ( task_name, meta['calls'][task_name][0]['executionStatus'], start, end )

    # create sorted table
    result = []
    for start in reversed(sorted(tasks.keys())):
        result.append(tasks[start])

    response["result"] = result
    if verbose:
        print(json.dumps(response, indent=4, sort_keys=True))
    return response

# simplifies the result JSON from that returned by cromwell
@jaws.command()
@click.argument('run_id')
def get_labels(run_id):
    """
    Get labels for a run.
    """
    url = config["CROMWELL"]['url'] + '/' + run_id + '/labels'
    r = requests.get(url)
    response = {}
    if r.status_code == 200:
        result = r.json()
        labels = result['labels']
        del labels['cromwell-workflow-id']
        if 'username' in labels:
            del labels['username']
        response["result"] = labels
    else:
        response["error"] = {
            "code" : 404,
            "message" : "Run not found"
        }
    if verbose:
        print(json.dumps(response, indent=4, sort_keys=True))
    return response

@jaws.command()
@click.option('--username', default=os.environ["USER"], help="Username; default=\$USER")
def queue(username):
    """
    Get a user's unfinished runs.
    """
    status = [ 'Submitted', 'Running', 'Aborting' ]
    labels = [ 'username:' + username ]
    data = {
        "status" : status,
        "label" : labels
    }
    url = config["CROMWELL"]["url"] + "/query"
    r = requests.get(url, data)
    response = {}
    if r.status_code == 200:
        response["result"] = r.json()
    else:
        response["error"] = {
            "code" : 404,
            "message" : "Run not found"
        }
    if verbose:
        print(json.dumps(response, indent=4, sort_keys=True))
    return response

@jaws.command()
@click.option('--username', default=os.environ["USER"], help="Username; default=\$USER")
@click.option('--days', default=10, help="Delta days, default=10")
def history(username, days):
    """
    Get a user's run history.
    """
    d = datetime.today() - timedelta(int(days))
    start_date = d.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    labels = [ "username:" + username ]
    data = {
        "label" : labels,
        "start" : start_date
    }
    url = config["CROMWELL"]['url'] + "/query"
    r = requests.get(url, data)
    response = {}
    if r.status_code == 200:
        response["result"] = r.json()
    else:
        response["error"] = {
            "code" : 404,
            "message" : "Run not found"
        }
    if verbose:
        print(json.dumps(response, indent=4, sort_keys=True))
    return response

@jaws.command()
@click.argument('json_file')
@click.argument('wdl_file')
@click.option('--zip_file', default=None, help="subworkflows zip")
@click.option('--username', default=os.environ["USER"], help="Username; default=\$USER")
@click.option('--labels', default=None, help="JSON file of labels")
def submit(wdl_file, json_file, username, labels, zip_file):
    """
    Submit a run.
    """
    tmp_basename = os.path.join(config["SITE"]["tmpdir"], str(uuid.uuid4()))
    if verbose:
        print("Using tmp base: %s" % (tmp_basename,)) # ECCE
        print("username : %s" % (username,)) # ECCE

    # labels
    tmp_labels_file = tmp_basename + ".labels.json"
    _prepare_labels_file(tmp_labels_file, username, labels)

    #_prepare_wdl_files()
    #_prepare_inputs_json()

    # SUBMIT TO CROMWELL
    files = {}
    files['workflowInputs'] = ( 'workflowInputs', open(json_file, 'r'), 'application/json' )
    files['workflowSource'] = ( 'workflowSource', open(wdl_file, 'r'), 'application/json' )
    files['labels'] = ( 'labels', open(tmp_labels_file, 'r'), 'application/json' )
    if zip_file is not None:
        files['workflowDependencies'] = ( 'workflowDependencies', open(zip_file, 'rb'), 'application/zip' )

    url = config["CROMWELL"]["url"]
    if verbose: print("Submitting run to %s" % (url,))
    r = requests.post(url, files=files)

    response = {}
    if r.status_code == 201:
        response["result"] = r.json()
    elif r.status_code == 500:
        response["error"] = {
            "code" : r.status_code,
            "message" : "Internal error"
        }
    elif r.status_code == 400:
        response["error"] = {
            "code" : r.status_code,
            "message" : "Invalid submission request"
        }
    else:
        response["error"] = {
            "code" : r.status_code,
            "message" : "Unknown error"
        }

    # CLEANUP
    os.remove(tmp_labels_file)

    if verbose:
        print(json.dumps(response, indent=4, sort_keys=True))
    return response

def _prepare_wdl_files(main_wdl_file, *wdl_files):
    """
    Validate WDL using womtool, identify File arguments in inputs JSON, remove any forbidden tags (e.g. "backend"), and create zip of subs.
    """ 
    sys.exit("NYI") # TODO

def _prepare_inputs_json(outfile, template, infile):
    """
    If inputs JSON file contains relative paths, make them absolute relative to the JSON file.  An inputs template JSON is required to identify which items are paths.
    """
    sys.exit("NYI") # TODO

def _prepare_labels_file(outfile, username, infile=None):
    """
    Given a username and optionally other labels (in a dict), create a JSON file.
    """
    assert outfile
    assert username
    labels = {}
    if infile is not None:
        with open(infile, 'r') as f:
            data = json.load(f)
            for key in data:
                labels[key] = data[key]
    labels["username"] = username
    with open(outfile, 'w') as f:
        json.dump(labels, f)

if __name__ == "__main__":
    """
    Main entry point
    """
    verbose=True # pretty-print JSON reponse
