"""
JAWS Analysis/Run management functions; these interact via REST with the JAWS Central server.
"""

import sys
import os
import json
import requests
import click
import logging
import uuid
import globus_sdk
from . import config, user, workflow


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def run():
    """JAWS Run-Workflows Commands"""
    pass


@run.command()
def queue():
    """List user's unfinished runs.

    :return: List of user's current runs in JSON format
    :rtype: str
    """
    current_user = user.User()
    url = f'{config.JawsConfig().get("JAWS", "url")}/run'
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@run.command()
@click.option("--days", default=1)
def history(days):
    """Print a list of the user's past runs.

    :param days: Time window to search, in days.
    :type days: int, optional
    """
    data = {"delta_days": days}
    url = f'{config.conf.get("JAWS", "url")}/search'
    current_user = user.User()
    try:
        r = requests.post(url, data=data, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


def _run_status(run_id):
    """Return the status of a run in JSON format.

    :param run_id: JAWS run ID
    :type run_id: int
    :return: Status of the run, in JSON string.
    :rtype: str
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}'
    current_user = user.User()
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    return result


@run.command()
@click.argument("run_id")
def status(run_id):
    """Print the current status of a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    result = _run_status(run_id)
    print(json.dumps(result, indent=4, sort_keys=True))


@run.command()
@click.argument("run_id")
def tasks(run_id):
    """Show status of each task of a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}'
    current_user = user.User()
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    metadata = r.json()
    status = metadata["status"]
    print("status: " + status)

    # failures
    if "failures" in metadata:
        print("failures:")
        for failure in metadata["failures"]:
            print("\t" + failure["message"])

    # call summary
    print("calls:")
    calls = metadata["calls"]
    result = []
    for task_name in calls.keys():
        task = calls[task_name]
        for attempt in task:
            status = attempt["executionStatus"]
            shard = attempt["shardIndex"]
            result.append((task_name, shard, status))
            if shard > -1:
                print(task_name + "[" + str(shard) + "]\t" + status)
            else:
                print("\t".join((task_name, status)))
            if status == "Failed":
                for failure in attempt["failures"]:
                    print("\t" + failure["message"])


@run.command()
@click.argument("run_id")
def metadata(run_id):
    """
    Print the detailed metadata for a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/metadata'
    current_user = user.User()
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@run.command()
@click.argument("run_id")
def log(run_id):
    """View the Cromwell log of a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/log'
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    print(r.text)


@run.command()
@click.argument("run_id")
def cancel(run_id):
    """Cancel a run; prints whether aborting was successful or not.

    :param run_id: JAWS run ID to cancel.
    :type run_id: int
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}'
    current_user = user.User()
    try:
        r = requests.delete(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code == 200:
        print("run " + run_id + " was canceled")
    else:
        sys.exit(r.text)


@run.command()
@click.argument("run_id")
@click.option("--task", default=None)
def delete(run_id, task):
    """Delete the output of a run or task to avoid caching.

    :param run_id: JAWS run ID
    :type run_id: int
    :param task: Name of the task to invalidate
    :type task: str, optional
    :return:
    """
    url = None
    if task is not None:
        url = f'{config.conf.get("JAWS", "url")}/run/{run_id}'
    else:
        url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/{task}'
    current_user = user.User()
    try:
        r = requests.delete(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    result = r.json()
    if r.status_code != 200:
        sys.exit(r.text)
    print(json.dumps(result, indent=4, sort_keys=True))


def _labels_string_to_dict(labels_string):
    """
    Converts string to dictionary.
    """
    labels_dict = {}
    labels = labels_string.split(",")
    for label in labels:
        if label.count(":") != 1:
            sys.exit(
                "Invalid label: " + label + '; labels must be in "key:value" format'
            )
        (key, value) = label.split(":")
        labels_dict[key] = value
    return labels_dict


# TODO add option for labels
@run.command()
@click.argument("wdl_file", nargs=1)
@click.argument("infile", nargs=1)
@click.argument("outdir", nargs=1)
@click.argument("site", nargs=1)
def submit(wdl_file, infile, outdir, site):
    """Submit a run for execution at a JAWS-Site.

    :param wdl_file: Path to workflow specification (WDL) file
    :type wdl_file: str
    :param infile: Path to inputs (JSON) file
    :type infile: str
    :param outdir: Path to output directory; doesn't have to exist
    :type outdir: str
    :param site: JAWS Site ID at which to run
    :type site: str
    """
    current_user = user.User()
    logger = logging.getLogger(__package__)

    # CONFIG
    conf = config.JawsConfig()
    JAWS_SITE = conf.get("JAWS", "site_id")
    GLOBUS_ENDPOINT = conf.get("GLOBUS", "endpoint")
    GLOBUS_BASEDIR = conf.get("GLOBUS", "basedir")
    JAWS_STAGING_SUBDIR = conf.get("JAWS", "staging_subdir")

    # DEFINE OUTPUT DIR AND VERIFY IT'S ACCESSIBLE VIA GLOBUS
    if not outdir:
        if current_user.config["USER"]["default_outdir"]:
            outdir = current_user.config["USER"]["default_outdir"]
        else:
            sys.exit(
                '--outdir required as no default specified in your config file' +
                '\nYou may set your "default_outdir" by editing %s'
                % (current_user.config_file,)
            )
    if not GLOBUS_BASEDIR:
        sys.exit("ERROR: $GLOBUS_BASEDIR env var not defined")
    if not outdir.startswith(GLOBUS_BASEDIR):
        sys.exit(
            "Invalid outdir, %s; Globus can only write under %s"
            % (outdir, GLOBUS_BASEDIR)
        )

    # CREATE UNIQUE STAGING ID
    submission_uuid = str(uuid.uuid4())

    # VALIDATE RUN
    run = workflow.Workflow(wdl_file, infile)
    if not run.validate():
        sys.exit("Failed validation; aborting")

    # PREPARE RUN
    staging_dir = os.path.join(GLOBUS_BASEDIR, JAWS_STAGING_SUBDIR)
    if not os.path.isdir(staging_dir):
        os.makedirs(staging_dir)
    run.prepare_wdls(staging_dir, submission_uuid)
    run.prepare_inputs(GLOBUS_BASEDIR, JAWS_STAGING_SUBDIR, JAWS_SITE, submission_uuid)

    # GET COMPUTE SITE INFO (E.G. GLOBUS ENDPOINT PARAMS)
    data = {"site": site, "max_ram_gb": run.max_ram_gb, "transfer_gb": run.transfer_gb}
    url = f'{config.conf.get("JAWS", "url")}/run/get_site'
    try:
        r = requests.post(url, data=data, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != requests.codes.ok:
        sys.exit(r.text)
    result = r.json()
    site = result["site"]
    print("Sending to: %s" % (site,))
    dest_endpoint = result["endpoint"]
    dest_dir = result["staging"]

    # GLOBUS TRANSFER
    current_user = user.User()
    transfer_client = current_user.transfer_client()
    tdata = globus_sdk.TransferData(
        transfer_client,
        GLOBUS_ENDPOINT,
        dest_endpoint,
        label=submission_uuid,
        sync_level="checksum",
        verify_checksum=True,
        preserve_timestamp=True,
        notify_on_succeeded=False,
        notify_on_failed=True,
        notify_on_inactive=True,
        skip_activation_check=False,
    )
    for source_file, dest_file in run.manifest:
        abs_dest_file = os.path.join(dest_dir, dest_file)
        if os.path.isdir(source_file):
            tdata.add_item(source_file, abs_dest_file, recursive=True)
        else:
            tdata.add_item(source_file, abs_dest_file, recursive=False)
    transfer_result = transfer_client.submit_transfer(tdata)
    transfer_task_id = transfer_result["task_id"]
    logger.info("Globus transfer task_id = %s" % (transfer_task_id,))

    # SUBMIT RUN
    # NOTE THAT THE FILE TRANSFER IS NOT COMPLETE YET
    data = {
        "site": site,
        "submission_uuid": submission_uuid,
        "transfer_task_id": transfer_task_id,
        "globus_endpoint": GLOBUS_ENDPOINT,
        "outdir": outdir,
    }
    url = f'{config.conf.get("JAWS", "url")}/run'
    logger.info("Submitting to %s:\n%s" % (url, data))
    try:
        r = requests.post(url, data=data, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != requests.codes.ok:
        sys.exit(r.text)
    result = r.json()
    run_id = result["run_id"]
    print("Successfully queued run %s" % (run_id,))
