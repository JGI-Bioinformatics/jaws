"""
JAWS Analysis/Run management functions; these interact via REST with the JAWS Central server.
"""

import os
import json
import requests
import click
import logging
import uuid
import shutil
from typing import Dict
from collections import defaultdict

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
        raise SystemExit(
            "Unable to communicate with JAWS server (too many redirects; bad url?)", err
        )
    except requests.exceptions.HTTPError as err:
        raise SystemExit("Unable to communicate with JAWS server (http error)", err)
    except requests.exceptions.RequestException as err:
        raise SystemExit("Unable to communicate with JAWS server", err)
    if r.status_code != 200:
        try:
            result = r.json()
        except Exception:
            raise SystemExit(r.text)
        if "detail" in result:
            raise SystemExit(result["detail"])
        else:
            raise SystemExit(result)
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
    url = f'{config.conf.get("JAWS", "url")}/search'
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
    if days < 1:
        raise SystemExit("User error: --days must be a positive integer")
    url = f'{config.conf.get("JAWS", "url")}/search/{days}'
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
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}'
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
@click.option("--fmt", default="text", help="Output format: text|json")
def task_status(run_id: int, fmt: str) -> None:
    """Show the current status of each task.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/task_status'
    r = _get(url)
    result = r.json()
    if fmt == "json":
        print(json.dumps(result, indent=4, sort_keys=True))
    else:
        print(
            "#TASK_NAME\tATTEMPT\tCROMWELL_JOB_ID\tSTATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON"
        )
        for log_entry in result:
            log_entry[1] = str(log_entry[1])
            log_entry[2] = str(log_entry[2])
            print("\t".join(log_entry))


@run.command()
@click.argument("run_id")
def metadata(run_id: int) -> None:
    """
    Print the detailed metadata for a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/metadata'
    r = _get(url)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@run.command()
@click.argument("run_id")
@click.option("--fmt", default="text", help="Output format: text|json")
def log(run_id: int, fmt: str) -> None:
    """View the log of Run state transitions.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/run_log'
    r = _get(url)
    result = r.json()
    if fmt == "json":
        print(json.dumps(result, indent=4, sort_keys=True))
    else:
        print("#STATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON")
        for log_entry in result:
            print("\t".join(log_entry))


@run.command()
@click.argument("run_id")
@click.option("--fmt", default="text", help="Output format: text|json")
def task_log(run_id: int, fmt: str) -> None:
    """Get log of each Tasks' state transitions.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    r = _get(f'{config.conf.get("JAWS", "url")}/run/{run_id}/task_log')
    result = r.json()
    if fmt == "json":
        print(json.dumps(result, indent=4, sort_keys=True))
    else:
        tasks = defaultdict(list)
        for log_entry in result:
            task_name = log_entry[0]
            tasks[task_name].append(log_entry)
        print(
            "#TASK_NAME\tATTEMPT\tCROMWELL_JOB_ID\tSTATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON"
        )
        for task_name in tasks:
            for log_entry in tasks[task_name]:
                log_entry[1] = str(log_entry[1])
                log_entry[2] = str(log_entry[2])
                print("\t".join(log_entry))


@run.command()
@click.argument("run_id")
@click.option("--failed", is_flag=True)
def output(run_id: int, failed: bool = False) -> None:
    """View the stdout/stderr output of Tasks.

    :param run_id: JAWS run ID
    :type run_id: int
    :param failed: Get output of failed tasks only
    :type failed: bool
    :return:
    """
    if failed:
        url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/failed'
    else:
        url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/output'
    r = _get(url)
    if r.status_code != 200:
        raise SystemExit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@run.command()
@click.argument("run_id")
def cancel(run_id):
    """Cancel a run; prints whether aborting was successful or not.

    :param run_id: JAWS run ID to cancel.
    :type run_id: int
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/cancel'
    current_user = user.User()
    try:
        r = requests.put(url, headers=current_user.header())
    except requests.exceptions.RequestException:
        raise SystemExit("Unable to communicate with JAWS server")
    result = r.json()
    if r.status_code != 201:
        if "detail" in result:
            raise SystemExit(result["detail"])
        else:
            raise SystemExit(r.text)
    print(json.dumps(result, indent=4, sort_keys=True))


def _list_sites() -> None:
    """List available Sites."""
    url = f'{config.conf.get("JAWS", "url")}/site'
    r = _get(url)
    result = r.json()
    print("Available Sites:")
    for a_site_id in result:
        print(f"  - {a_site_id}")


@run.command()
def list_sites() -> None:
    """List available Sites"""
    _list_sites()


@run.command()
@click.argument("wdl_file", nargs=1)
@click.argument("infile", nargs=1)
@click.argument("outdir", nargs=1)
@click.argument("site", nargs=1)
@click.option("--out_endpoint", default=None, help="Globus endpoint to send output")
def submit(wdl_file, infile, outdir, site, out_endpoint):
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
    logger = logging.getLogger(__package__)

    current_user = user.User()

    # STAGING DIR
    staging_subdir = config.Configuration().get("USER", "staging_dir")
    globus_basedir = config.Configuration().get("GLOBUS", "basedir")
    if not staging_subdir.startswith(globus_basedir):
        raise SystemExit(
            f"Configuration error: Staging dir must be under endpoint's basedir: {globus_basedir}"
        )
    if not os.path.isdir(staging_subdir):
        os.makedirs(staging_subdir)

    # GLOBUS
    local_endpoint_id = config.conf.get("GLOBUS", "endpoint_id")
    if out_endpoint is None:
        out_endpoint = local_endpoint_id

    # OUTDIR
    outdir = os.path.abspath(outdir)
    if out_endpoint == local_endpoint_id and not outdir.startswith(globus_basedir):
        raise SystemExit(f"Outdir must be under endpoint's basedir: {globus_basedir}")
    if not os.path.isdir(outdir):
        try:
            os.makedirs(outdir)
        except Exception as error:
            raise SystemExit(f"Unable to make outdir, {outdir}: {error}")

    # CONFIRM OUTDIR WRITABLE BY COPYING WDL, INPUT FILES
    try:
        shutil.copy2(wdl_file, outdir)
    except Exception as error:
        raise SystemExit(f"Unable to copy wdl file to outdir: {error}")
    try:
        shutil.copy2(infile, outdir)
    except Exception as error:
        raise SystemExit(f"Unable to copy json file to outdir: {error}")

    # GET SITE INFO
    compute_site_id = site.upper()
    url = f'{config.conf.get("JAWS", "url")}/site/{compute_site_id}'
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.exceptions.RequestException:
        raise SystemExit("Unable to communicate with JAWS server")
    if r.status_code == 404:
        print(f"{compute_site_id} is not a valid Site ID.")
        _list_sites()
        raise SystemExit("Please try again with a valid Site ID")
    elif r.status_code != requests.codes.ok:
        result = r.json()
        raise SystemExit(result["detail"])
    result = r.json()
    compute_basedir = result["globus_basepath"]
    compute_staging_subdir = result["staging_subdir"]
    compute_max_ram_gb = int(result["max_ram_gb"])

    # VALIDATE WORKFLOW
    submission_id = str(uuid.uuid4())
    wdl = workflow.WdlFile(wdl_file, submission_id)
    inputs_json = workflow.WorkflowInputs(infile, submission_id)

    jaws_site_staging_dir = workflow.join_path(compute_basedir, compute_staging_subdir)
    local_staging_endpoint = workflow.join_path(globus_basedir, staging_subdir)
    manifest_file = workflow.Manifest(local_staging_endpoint, jaws_site_staging_dir)

    wdl.validate()
    inputs_json.validate()

    max_ram_gb = wdl.max_ram_gb
    if max_ram_gb > compute_max_ram_gb:
        raise AnalysisError(
            f"The workflow requires {max_ram_gb}GB but {compute_site_id} has only {compute_max_ram_gb}GB available"
        )

    site_id = config.conf.get("JAWS", "site_id")
    site_subdir = workflow.join_path(local_staging_endpoint, site_id)

    sanitized_wdl, zip_file = workflow.compress_wdls(wdl, local_staging_endpoint)
    moved_files = workflow.move_input_files(inputs_json, site_subdir)

    staged_json = workflow.join_path(local_staging_endpoint, f"{submission_id}.json")
    jaws_site_staging_site_subdir = workflow.join_path(jaws_site_staging_dir, site_id)
    modified_json = inputs_json.prepend_paths_to_json(jaws_site_staging_site_subdir)
    modified_json.write_to(staged_json)

    manifest_file.add(sanitized_wdl, zip_file, staged_json, *moved_files)
    staged_manifest = workflow.join_path(staging_subdir, f"{submission_id}.tsv")
    manifest_file.write_to(staged_manifest)

    # SUBMIT RUN TO CENTRAL
    data = {
        "site_id": compute_site_id,
        "submission_id": submission_id,
        "input_site_id": config.conf.get("JAWS", "site_id"),
        "input_endpoint": config.conf.get("GLOBUS", "endpoint_id"),
        "output_endpoint": out_endpoint,
        "output_dir": outdir,
    }
    files = {"manifest": open(staged_manifest, "r")}
    url = f'{config.conf.get("JAWS", "url")}/run'
    logger.debug(f"Submitting run: {data}")
    try:
        r = requests.post(url, data=data, files=files, headers=current_user.header())
    except requests.exceptions.RequestException:
        raise SystemExit("Unable to communicate with JAWS server")
    result = r.json()
    if r.status_code != 201:
        raise SystemExit(result["detail"])
    if "run_id" not in result:
        raise SystemExit(f"Run submission failed: {result}")
    run_id = result["run_id"]
    logger.info(f"Submitted run {run_id}: {data}")

    # WRITE RUN INFO TO OUTDIR
    logfile = os.path.join(outdir, f"jaws.run_{run_id}.log")
    with open(logfile, "w") as fh:
        fh.write(r.text + "\n")
    print(r.text)


@run.command()
@click.argument("wdl_file", nargs=1)
def inputs(wdl_file: str) -> None:
    """Generate inputs template (JSON) from workflow (WDL) file.

    :param wdl_file: Path to workflow specification (WDL) file
    :type wdl_file: str
    :return:
    """
    if not os.path.isfile(wdl_file):
        raise IOError(f"File not found: {wdl_file}")
    stdout, stderr = workflow.womtool("inputs", wdl_file)
    if stderr:
        raise SystemExit(stderr)
    print(stdout.strip())


@run.command()
@click.argument("wdl_file", nargs=1)
def validate(wdl_file: str) -> None:
    """Validate a WDL using Cromwell's WOMTool.

    :param wdl_file: Path to workflow specification (WDL) file
    :type wdl_file: str
    :return:
    """
    if not os.path.isfile(wdl_file):
        raise IOError(f"File not found: {wdl_file}")
    stdout, stderr = workflow.womtool("inputs", wdl_file)
    if stderr:
        raise SystemExit(stderr)
    else:
        print("Workflow is OK")
