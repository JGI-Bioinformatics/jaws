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
            "#CROMWELL_RUN_ID\tTASK_NAME\tATTEMPT\tCROMWELL_JOB_ID\tSTATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON"
        )
        for row in result:
            row[2] = str(row[2])
            row[3] = str(row[3])
            print("\t".join(row))


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
        print(
            "#CROMWELL_RUN_ID\tTASK_NAME\tATTEMPT\tCROMWELL_JOB_ID\tSTATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON"
        )
        for row in result:
            row[2] = str(row[2])
            row[3] = str(row[3])
            print("\t".join(row))


@run.command()
@click.argument("run_id")
@click.option("--fmt", default="text", help="Output format: text|json")
def errors(run_id: int, fmt: str) -> None:
    """View error messages and stderr for failed tasks.

    :param run_id: JAWS run ID
    :type run_id: int
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/errors'
    r = _get(url)
    if r.status_code != 200:
        raise SystemExit(r.text)
    result = r.json()
    if fmt == "json":
        print(json.dumps(result, indent=4, sort_keys=True))
    else:
        for task_name in result:
            print(f"{task_name}:")
            print(result[task_name])
            print("\n")


@run.command()
@click.argument("run_id")
def cancel(run_id):
    """Cancel a particular run.

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


@run.command()
def cancel_all():
    """Cancel all active runs."""
    url = f'{config.conf.get("JAWS", "url")}/run/cancel-all'
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
@click.argument("json_file", nargs=1)
@click.argument("site", nargs=1)
@click.option("--tag", default="")
def submit(wdl_file: str, json_file: str, site: str, tag: str):
    """Submit a run for execution at a JAWS-Site.

    :param wdl_file: Path to workflow specification (WDL) file
    :type wdl_file: str
    :param json_file: Path to input JSON file
    :type json_file: str
    :param site: JAWS Site ID at which to run
    :type site: str
    :param tag: User-supplied label for this run.
    :type tag: str
    """
    logger = logging.getLogger(__package__)

    # the users' jaws id may not match the linux uid where the client is installed
    current_user = user.User()
    user_url = f'{config.conf.get("JAWS", "url")}/user'
    user_rec = _get(user_url)
    user_json = user_rec.json()
    uid = user_json["uid"]

    staging_subdir = config.Configuration().get("JAWS", "staging_dir")
    staging_user_subdir = os.path.join(staging_subdir, uid)
    globus_host_path = config.Configuration().get("GLOBUS", "host_path")
    output_directory = config.conf.get("JAWS", "data_repo_basedir")
    input_site_id = config.conf.get("JAWS", "site_id")
    local_staging_endpoint = workflow.join_path(globus_host_path, staging_user_subdir)

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
    compute_basedir = result["globus_host_path"]
    compute_uploads_subdir = result["uploads_dir"]
    compute_max_ram_gb = int(result["max_ram_gb"])

    # VALIDATE WORKFLOW WDLs
    submission_id = str(uuid.uuid4())
    try:
        wdl = workflow.WdlFile(wdl_file, submission_id)
    except workflow.WdlError as error:
        raise SystemExit(f"There is a problem with your workflow:\n{error}")
    try:
        wdl.validate()
    except workflow.WdlError as error:
        raise SystemExit(error)
    max_ram_gb = wdl.max_ram_gb
    if max_ram_gb > compute_max_ram_gb:
        raise SystemExit(
            f"The workflow requires {max_ram_gb}GB but {compute_site_id} has only {compute_max_ram_gb}GB available"
        )

    # any and all subworkflow WDL files must be supplied to Cromwell in a single ZIP archive
    try:
        staged_wdl, zip_file = workflow.compress_wdls(wdl, local_staging_endpoint)
    except IOError as error:
        raise SystemExit(f"Unable to copy WDLs to inputs dir: {error}")

    # VALIDATE INPUTS JSON
    try:
        inputs_json = workflow.WorkflowInputs(json_file, submission_id)
    except json.JSONDecodeError as error:
        raise SystemExit(f"Your file, {json_file}, is not a valid JSON file: {error}")

    staged_json = workflow.join_path(local_staging_endpoint, f"{submission_id}.json")
    site_subdir = workflow.join_path(local_staging_endpoint, input_site_id)
    jaws_site_staging_dir = workflow.join_path(compute_basedir, compute_uploads_subdir)
    jaws_site_staging_site_subdir = workflow.join_path(
        jaws_site_staging_dir, input_site_id
    )

    # copy infiles in inputs json to site's inputs dir so they may be read by jaws user and
    # transferred to the compute site via globus
    moved_files = inputs_json.move_input_files(site_subdir)

    # the paths in the inputs json file are changed to their paths at the compute site
    modified_json = inputs_json.prepend_paths_to_json(jaws_site_staging_site_subdir)
    modified_json.write_to(staged_json)

    # the original inputs json file is kept with the run submission for record-keeping only
    orig_json = workflow.join_path(local_staging_endpoint, f"{submission_id}.orig.json")
    try:
        shutil.copy(json_file, orig_json)
    except IOError as error:
        raise SystemExit(f"Error copying JSON to {orig_json}: {error}")
    try:
        os.chmod(orig_json, 0o0664)
    except PermissionError as error:
        raise SystemExit(f"Unable to chmod {orig_json}: {error}")

    # write the file transfer manifest; jaws-central shall submit the transfer to globus
    manifest_file = workflow.Manifest(local_staging_endpoint, compute_uploads_subdir)
    manifest_file.add(staged_wdl, zip_file, staged_json, orig_json, *moved_files)
    staged_manifest = workflow.join_path(staging_user_subdir, f"{submission_id}.tsv")
    manifest_file.write_to(staged_manifest)

    # SUBMIT RUN TO CENTRAL
    local_endpoint_id = config.conf.get("GLOBUS", "endpoint_id")
    data = {
        "site_id": compute_site_id,
        "submission_id": submission_id,
        "input_site_id": input_site_id,
        "input_endpoint": local_endpoint_id,
        "output_endpoint": local_endpoint_id,  # return to original submission site
        "output_dir": output_directory,  # jaws-writable dir to initially receive results
        "wdl_file": wdl_file,
        "json_file": json_file,
        "tag": tag,
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


@run.command()
@click.argument("run_id")
@click.argument("dest")
def get(run_id: int, dest: str) -> None:
    """Copy the output of a run to the specified folder.

    :param run_id: JAWS run ID
    :type run_id: int
    :param dest: destination path
    :type dest: str
    :return:
    """
    logger = logging.getLogger(__package__)
    result = _run_status(run_id)
    status = result["status"]
    src = result["output_dir"]

    if status != "download complete":
        raise SystemExit(
            f"Run {run_id} output is not yet available; status is {status}"
        )

    if src is None:
        logger.error(f"Run {run_id} doesn't have an output_dir defined")
        raise SystemExit(f"Run {run_id} doesn't have an output_dir defined")

    try:
        result = workflow.rsync(
            src,
            dest,
            ["-rLtq", "--chmod=Du=rwx,Dg=rwx,Do=,Fu=rw,Fg=rw,Fo=", "--exclude=inputs"],
        )
    except IOError as error:
        logger.error(f"Rsync output failed for run {run_id}: {error}")
        raise SystemExit(f"Error getting output for run {run_id}: {error}")
    if result.returncode != 0:
        err_msg = f"Failed to rsync {src}->{dest}: {result.stdout}; {result.stderr}"
        raise SystemExit(err_msg)
