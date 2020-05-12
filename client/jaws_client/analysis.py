"""
JAWS Analysis/Run management functions; these interact via REST with the JAWS Central server.
"""

import os
import json
import requests
import click
import logging
import subprocess
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
        raise SystemExit(r.text)
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
def tasks(run_id: int) -> None:
    """Show status of each task of a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/tasks'
    r = _get(url)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


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
def logs(run_id: int) -> None:
    """View the Cromwell logs of a run.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/logs'
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
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}/abort'
    current_user = user.User()
    try:
        r = requests.put(url, headers=current_user.header())
    except requests.exceptions.RequestException:
        raise SystemExit("Unable to communicate with JAWS server")
    if r.status_code != 201:
        raise SystemExit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@run.command()
@click.argument("run_id")
def uncache(run_id: int) -> None:
    """Prevent the results of a run from being reused by the cache.

    :param run_id: JAWS run ID
    :type run_id: int
    :return:
    """
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/run/{run_id}'
    try:
        r = requests.delete(url, headers=current_user.header())
    except requests.exceptions.RequestException:
        raise SystemExit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        raise SystemExit(r.text)
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
    outdir = os.path.abspath(outdir)

    compute_site_id = site.upper()
    current_user = user.User()
    logger = logging.getLogger(__package__)

    globus_basedir = config.Configuration().get("GLOBUS", "basedir")
    staging_subdir = config.Configuration().get("USER", "staging_dir")

    if not staging_subdir.startswith(globus_basedir):
        raise SystemExit(
            f"Staging dir must be under endpoint's basedir: {globus_basedir}"
        )
    if not os.path.isdir(staging_subdir):
        os.makedirs(staging_subdir)
    local_endpoint_id = config.conf.get("GLOBUS", "endpoint_id")

    if out_endpoint is None:
        out_endpoint = local_endpoint_id
    if out_endpoint == local_endpoint_id and not outdir.startswith(globus_basedir):
        raise SystemExit(f"Outdir must be under endpoint's basedir: {globus_basedir}")
    if not os.path.isdir(outdir):
        try:
            os.makedirs(outdir)
        except Exception as error:
            raise SystemExit(f"Invalid outdir, {outdir}: {error}")
    try:
        shutil.copy2(wdl_file, outdir)
    except Exception as error:
        raise SystemExit(f"Unable to copy wdl file to outdir: {error}")
    try:
        shutil.copy2(infile, outdir)
    except Exception as error:
        raise SystemExit(f"Unable to copy json file to outdir: {error}")

    submission_id = str(uuid.uuid4())

    url = f'{config.conf.get("JAWS", "url")}/site/{compute_site_id}'

    try:
        r = requests.get(url, headers=current_user.header())
    except requests.exceptions.RequestException:
        raise SystemExit("Unable to communicate with JAWS server")
    if r.status_code == 404:
        print(f"{compute_site_id} is not a valid Site ID.")
        list_sites()
        raise SystemExit("Please try again with a valid Site ID")
    elif r.status_code != requests.codes.ok:
        raise SystemExit(r.text)

    result = r.json()
    compute_basedir = result["globus_basepath"]
    compute_staging_subdir = result["staging_subdir"]
    compute_max_ram_gb = int(result["max_ram_gb"])

    # PREPARE RUN
    jaws_site_staging_dir = workflow.join_path(compute_basedir, compute_staging_subdir)
    local_staging_endpoint = workflow.join_path(globus_basedir, staging_subdir)

    wdl = workflow.WdlFile(wdl_file, submission_id)
    inputs_json = workflow.WorkflowInputs(infile, submission_id)
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

    compressed_wdl, zip_file = workflow.compress_wdls(wdl, local_staging_endpoint)
    moved_files = workflow.move_input_files(inputs_json, site_subdir)

    staged_json = workflow.join_path(local_staging_endpoint, f"{submission_id}.json")
    modified_json = inputs_json.prepend_paths_to_json(jaws_site_staging_dir)
    modified_json.write_to(staged_json)

    manifest_file.add(compressed_wdl, zip_file, staged_json, *moved_files)
    staged_manifest = workflow.join_path(staging_subdir, f"{submission_id}.tsv")
    manifest_file.write_to(staged_manifest)

    data = {
        "site_id": site,
        "submission_id": submission_id,
        "input_site_id": config.conf.get("JAWS", "site_id"),
        "input_endpoint": config.conf.get("GLOBUS", "endpoint_id"),
        "output_endpoint": out_endpoint,
        "output_dir": outdir,
    }
    files = {"manifest": open(staged_manifest, "r")}
    url = f'{config.conf.get("JAWS", "url")}/run'
    logger.debug("Submitting run to %s:\n%s" % (url, data))
    current_user = user.User()
    try:
        r = requests.post(url, data=data, files=files, headers=current_user.header())
    except requests.exceptions.RequestException:
        raise SystemExit("Unable to communicate with JAWS server")
    if r.status_code != 201:
        raise SystemExit(r.text)
    result = r.json()
    if "run_id" not in result:
        raise SystemExit(f"Invalid response from JAWS: {result}")
    run_id = result["run_id"]
    logfile = os.path.join(outdir, f"jaws.run_{run_id}.log")
    with open(logfile, "w") as fh:
        fh.write(r.text+"\n")
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
    proc = subprocess.run(
        ["java", "-jar", config.conf.get("JAWS", "womtool_jar"), "inputs", wdl_file, ],
        capture_output=True,
        text=True,
    )
    if proc.stderr:
        raise SystemExit(proc.stderr)
    print(proc.stdout.strip())
