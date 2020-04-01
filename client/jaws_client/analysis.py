"""
JAWS Analysis/Run management functions; these interact via REST with the JAWS Central server.
"""

import sys
import os
import json
import requests
import click
import logging
from typing import Dict
from jaws_client import config, user, workflow


class AnalysisError(Exception):
    def __init__(self, message):
        super().__init__(message)


def _get(url):
    """REST GET.  Adds current user token to header.  Dies on error.

    :param url: URL
    :type url: str
    :return: requests object
    :rtype: requests
    """
    current_user = user.User()
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.exceptions.Timeout as err:
        raise SystemExit("Unable to communicate with JAWS server (timeout)", err)
    except requests.exceptions.TooManyRedirects as err:
        raise SystemExit("Unable to communicate with JAWS server (too many redirects; bad url?)", err)
    except requests.exceptions.HTTPError as err:
        raise SystemExit("Unable to communicate with JAWS server (http error)", err)
    except requests.exceptions.RequestException as err:
        raise SystemExit("Unable to communicate with JAWS server", err)
    if r.status_code != 200:
        sys.exit(r.text)
    return r


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def run():
    """JAWS Run-Workflows Commands"""
    pass


@run.command()
def queue() -> None:
    """List user's unfinished runs.

    :return: List of user's current runs in JSON format
    :rtype: str
    """
    url = f'{config.JawsConfig().get("JAWS", "url")}/search'
    r = _get(url)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@run.command()
@click.option("--days", default=1)
def history(days: int) -> None:
    """Print a list of the user's past runs.

    :param days: Time window to search, in days.
    :type days: int, optional
    """
    url = f'{config.JawsConfig().get("JAWS", "url")}/search/{days}'
    r = _get(url)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


def _run_status(run_id: int) -> Dict[str, str]:
    """Return the status of a run in JSON format.

    :param run_id: JAWS run ID
    :type run_id: int
    :return: Status of the run, in JSON string.
    :rtype: str
    """
    url = f'{config.JawsConfig().get("JAWS", "url")}/run/{run_id}'
    r = _get(url)
    result = r.json()
    return result


@run.command()
@click.argument("run_id")
def status(run_id: int) -> None:
    """Print the current status of a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    result = _run_status(run_id)
    print(json.dumps(result, indent=4, sort_keys=True))


@run.command()
@click.argument("run_id")
def tasks(run_id: int) -> None:
    """Show status of each task of a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.JawsConfig().get("JAWS", "url")}/run/{run_id}'
    r = _get(url)
    result = r.json()
    metadata = r.json()
    if "run_id" not in metadata:
        sys.exit(f"Invalid response from JAWS: {metadata}")
    status = metadata["status"]
    print("status: " + status)

    # failures
    if "failures" in metadata:
        print("failures:")
        for failure in metadata["failures"]:
            print("\t" + failure["message"])

    # call summary
    if "calls" not in metadata:
        return
    print("calls:")
    calls = metadata["calls"]
    result = []
    for task_name in calls.keys():
        task = calls[task_name]
        for attempt in task:
            if "executionStatus" not in attempt:
                continue
            status = attempt["executionStatus"]
            if "shardIndex" not in attempt:
                continue
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
def metadata(run_id: int) -> None:
    """
    Print the detailed metadata for a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.JawsConfig().get("JAWS", "url")}/run/{run_id}/metadata'
    r = _get(url)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@run.command()
@click.argument("run_id")
def logs(run_id: int) -> None:
    """View the Cromwell logs of a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.JawsConfig().get("JAWS", "url")}/run/{run_id}/logs'
    r = _get(url)
    print(r.text)


@run.command()
@click.argument("run_id")
def errors(run_id):
    """View the logs for failed tasks.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/errors'
    r = _get(url)
    print(r.text)


@run.command()
@click.argument("run_id")
def cancel(run_id):
    """Cancel a run; prints whether aborting was successful or not.

    :param run_id: JAWS run ID to cancel.
    :type run_id: int
    """
    url = f'{config.JawsConfig().get("JAWS", "url")}/run/{run_id}'
    _get(url)


@run.command()
@click.argument("run_id")
@click.option("--task", default=None)
def delete(run_id: int, task: str) -> None:
    """Delete the output of a run or task to avoid caching.

    :param run_id: JAWS run ID
    :type run_id: int
    :param task: Name of the task to invalidate
    :type task: str, optional
    :return:
    """
    url = None
    if task is not None:
        url = f'{config.JawsConfig().get("JAWS", "url")}/run/{run_id}'
    else:
        url = f'{config.JawsConfig().get("JAWS", "url")}/run/{run_id}/{task}'
    r = _get(url)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@run.command()
def list_sites() -> None:
    """List available Sites"""
    url = f'{config.conf.get("JAWS", "url")}/site'
    r = _get(url)
    result = r.json()
    print("Available Sites:")
    for a_site_id in result:
        print(f"  - {a_site_id}")


@run.command()
@click.argument("wdl_file", nargs=1)
@click.argument("infile", nargs=1)
@click.argument("outdir", nargs=1)
@click.argument("site", nargs=1)
@click.option("--out_ep", default=None, help="Globus endpoint to send output")
def submit(wdl_file: str, infile: str, outdir: str, site: str, out_ep: str) -> None:
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
    compute_site_id = site.upper()
    current_user = user.User()
    logger = logging.getLogger(__package__)

    globus_basedir = config.conf.get("GLOBUS", "basedir")
    staging_subdir = config.conf.get("JAWS", "staging_subdir")
    staging_dir = os.path.join(globus_basedir, staging_subdir)
    if not os.path.isdir(staging_dir):
        os.makedirs(staging_dir)

    local_endpoint_id = config.conf.get("GLOBUS", "endpoint_id")
    if out_ep is None:
        out_ep = local_endpoint_id
    if out_ep == local_endpoint_id and not outdir.startswith(globus_basedir):
        sys.exit(f"Outdir must be under endpoint's basedir: {globus_basedir}")

    # GET COMPUTE SITE INFO
    url = f'{config.conf.get("JAWS", "url")}/site/{compute_site_id}'
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.exceptions.RequestException:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code == 404:
        print(f"{compute_site_id} is not a valid Site ID.")
        list_sites()
        sys.exit("Please try again with a valid Site ID")
    elif r.status_code != requests.codes.ok:
        sys.exit(r.text)
    result = r.json()
    compute_basedir = result["globus_basepath"]
    compute_staging_subdir = result["staging_subdir"]
    compute_max_ram_gb = int(result["max_ram_gb"])

    # PREPARE RUN
    run = workflow.Workflow(wdl_file, infile)
    run.validate()
    max_ram_gb = run.max_ram_gb()
    if max_ram_gb > compute_max_ram_gb:
        raise AnalysisError(
            f"The workflow requires {max_ram_gb}GB but {compute_site_id} has only {compute_max_ram_gb}GB available"
        )
    submission_id = run.prepare_submission(compute_basedir, compute_staging_subdir)

    # SUBMIT RUN
    data = {
        "site_id": site,
        "submission_id": submission_id,
        "input_site_id": config.conf.get("JAWS", "site_id"),
        "input_endpoint": config.conf.get("GLOBUS", "endpoint_id"),
        "output_endpoint": out_ep,
        "output_dir": outdir,
    }
    files = {"manifest": open(run.manifest_file, "r")}
    url = f'{config.conf.get("JAWS", "url")}/run'
    logger.debug("Submitting run to %s:\n%s" % (url, data))
    current_user = user.User()
    try:
        r = requests.post(url, data=data, files=files, headers=current_user.header())
    except requests.exceptions.RequestException:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != requests.codes.ok:
        sys.exit(r.text)
    result = r.json()
    if "run_id" not in result:
        sys.exit(f"Invalid response from JAWS: {result}")
    run_id = result["run_id"]
    print(f"Successfully queued run {run_id}")
