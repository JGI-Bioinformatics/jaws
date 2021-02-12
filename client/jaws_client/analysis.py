"""
JAWS Analysis/Run management functions; these interact via REST with the JAWS Central server.
"""

import os
import json
import requests
import click
import logging
import uuid
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
            "#TASK_NAME\tATTEMPT\tCROMWELL_JOB_ID\tSTATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON\tSTATUS_DETAIL"
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
@click.argument("compute_site_id", nargs=1)
@click.option("--outep", default=None, help="Globus endpoint to send output")
@click.option("--outdir", default=None, help="Directory to send output")
def submit(wdl_file, infile, compute_site_id, outep, outdir):
    """Submit a run for execution at a JAWS-Site.

    :param wdl_file: Path to workflow specification (WDL) file
    :type wdl_file: str
    :param infile: Path to inputgetpass.getuser()s (JSON) file
    :type infile: str
    :param compute_site_id: JAWS Site ID at which to run
    :type compute_site_id: str
    :param output_endpoint: Globus endpoint id
    :type output_endpoint: str
    :param output_dir: Output dir, as expected by Globus
    :type output_dir: str
    """
    # short variable names were used for user convenience
    output_globus_endpoint = outep
    output_dir = outdir

    # The "output_dir" should be in the form expected by Globus.
    # Specifically, if the path starts with a '/' (looks like absolute path),
    # Globus shall prepend the Globus endpoint's host_path; otherwise
    # Globus shall prepend the Globus endpoint's default_dir.
    #
    # The Globus endpoint's host_path and default_dir are part of the endpoint's config,
    # and are neither provided nor needed by JAWS.
    #
    # For example, if endpoint has host_path='/scratch' and default_dir='~' then:
    # ex1. /mydir => /scratch/mydir
    # ex2. mydir => ~/mydir

    logger = logging.getLogger(__package__)

    current_user = user.User()

    # The "submission_id" is used to name folders.
    # Once the run has been submitted (successfully), an integer primary key,
    # "run_id", will be issued and returned to the user.
    # The user uses only the run_id to refer to the run.
    # That is, submission_id is only used internally.
    submission_id = str(uuid.uuid4())

    # GET INPUT (LOCAL) SITE INFO
    input_site_id = config.conf.get("JAWS", "site_id")
    url = f'{config.conf.get("JAWS", "url")}/site/{input_site_id}'
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.exceptions.RequestException:
        raise SystemExit("Unable to communicate with JAWS server")
    if r.status_code == 404:
        print(f"{input_site_id} is not a valid Site ID.")
        _list_sites()
        raise SystemExit("Please try again with a valid input Site ID")
    elif r.status_code != requests.codes.ok:
        result = r.json()
        raise SystemExit(result["detail"])
    result = r.json()
    input_globus_endpoint = result["globus_endpoint"]
    #input_globus_host_path = result["globus_host_path"]
    input_dir = result["input_dir"]

    # GET OUTPUT INFO (by default, also local)
    if output_globus_endpoint is None:
        # return results back to this Site
        output_globus_endpoint = input_globus_endpoint
        output_dir = os.path.join(result["output_dir"], submission_id)
    elif output_dir is None:
        # put in the endpoint's default dir
        output_dir = submission_id

    # GET COMPUTE SITE INFO
    compute_site_id = compute_site_id.upper()
    url = f'{config.conf.get("JAWS", "url")}/site/{compute_site_id}'
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.exceptions.RequestException:
        raise SystemExit("Unable to communicate with JAWS server")
    if r.status_code == 404:
        print(f"{compute_site_id} is not a valid Site ID.")
        _list_sites()
        raise SystemExit("Please try again with a valid compute Site ID")
    elif r.status_code != requests.codes.ok:
        result = r.json()
        raise SystemExit(result["detail"])
    result = r.json()
    #compute_host_path = result["globus_host_path"]
    compute_input_dir = result["input_dir"]
    compute_max_ram_gb = int(result["max_ram_gb"])

    # VALIDATE WORKFLOW AND INPUTS JSON
    try:
        wdl = workflow.WdlFile(wdl_file)
    except workflow.WdlError as error:
        raise SystemExit(f"There is a problem with your workflow:\n{error}")
    except Exception as error:
        raise SystemExit(f"Unexpected error validating workflow: {error}")
    try:
        inputs_json = workflow.WorkflowInputs(infile, submission_id)
    except Exception as error:
        raise SystemExit(f"Your file, {infile}, is not a valid JSON file: {error}")
    # check if valid wdl(s)
    wdl.validate()
    # check if valid json and all infiles accessible
    inputs_json.validate()
    # check if requested compute site has sufficient ram
    max_ram_gb = wdl.max_ram_gb
    if max_ram_gb > compute_max_ram_gb:
        raise AnalysisError(
            f"The workflow requires {max_ram_gb}GB but {compute_site_id} has only {compute_max_ram_gb}GB available"
        )

    # PREPARE ALL INPUT FILES

    # sanitize wdl and zip subworkflows
    sanitized_wdl, zip_file = workflow.compress_wdls(wdl, submission_id, input_dir)

    # copy input files to input dir, in subdirectory named after the input-Site
    input_files_dir = workflow.join_path(input_dir, input_site_id)
    staged_files = workflow.copy_input_files(inputs_json, input_files_dir)

    # write a new json file for the compute site; all the input file paths therein shall
    # match the locations of the destination files after having been transferred to
    # the compute-Site, in a subdirectory named after the input-Site
    compute_input_files_dir = os.path.join(compute_input_dir, input_site_id)
    compute_json = inputs_json.prepend_paths_to_json(compute_input_files_dir)
    compute_json_file = workflow.join_path(input_dir, f"{submission_id}.json")
    compute_json.write_to(compute_json_file)

    # write the Globus transfer manifest file, which is a table with source (input site) and
    # destination (compute site) paths
    manifest = workflow.Manifest(input_dir, compute_input_dir)
    manifest.add(sanitized_wdl, zip_file, compute_json_file, *staged_files)
    manifest_file = workflow.join_path(input_dir, f"{submission_id}.tsv")
    manifest.write_to(manifest_file)

    # SUBMIT RUN TO CENTRAL
    data = {
        "submission_id": submission_id,
        "input_site_id": input_site_id,
        "input_globus_endpoint": input_globus_endpoint,
        "compute_site_id": compute_site_id,
        "output_globus_endpoint": output_globus_endpoint,
        "output_dir": output_dir,
    }
    files = {"manifest": open(manifest_file, "r")}
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
