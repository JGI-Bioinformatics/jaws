"""
JAWS CLI
"""

import click
import os
import requests
import json
import uuid
import shutil
from typing import Dict
from jaws_client import log as logging
from jaws_client import deprecated
from jaws_client import wfcopy as wfc
from jaws_client.config import Configuration

JAWS_LOG_ENV = "JAWS_CLIENT_LOG"
JAWS_USER_LOG = os.path.expanduser("~/jaws.log")
JAWS_CONFIG_ENV = "JAWS_CLIENT_CONFIG"
JAWS_CONFIG_DEFAULT_FILE = os.path.expanduser("~/jaws-client.conf")
JAWS_USER_CONFIG_ENV = "JAWS_USER_CONFIG"
JAWS_USER_CONFIG_DEFAULT_FILE = os.path.expanduser("~/jaws.conf")


logger = None
config = None


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "jaws_config_file", default=None, help="JAWS config file")
@click.option("--user", "user_config_file", default=None, help="User config file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option(
    "--log-level",
    "log_level",
    default="INFO",
    help="Logging level [debug|info|warning|error|critical]",
)
def main(jaws_config_file: str, user_config_file: str, log_file: str, log_level: str):
    """JGI Analysis Workflows Service"""
    if log_file is None:
        log_file = (
            os.environ[JAWS_LOG_ENV] if JAWS_LOG_ENV in os.environ else JAWS_USER_LOG
        )
    global logger
    logger = logging.setup_logger(__package__, log_file, log_level)

    if jaws_config_file is None:
        jaws_config_file = (
            os.environ[JAWS_CONFIG_ENV]
            if JAWS_CONFIG_ENV in os.environ
            else JAWS_CONFIG_DEFAULT_FILE
        )
    os.environ[JAWS_CONFIG_ENV] = jaws_config_file
    if user_config_file is None:
        user_config_file = (
            os.environ[JAWS_USER_CONFIG_ENV]
            if JAWS_USER_CONFIG_ENV in os.environ
            else JAWS_USER_CONFIG_DEFAULT_FILE
        )
    os.environ[JAWS_USER_CONFIG_ENV] = user_config_file
    global config
    config = Configuration(jaws_config_file, user_config_file)


@main.command()
def health() -> None:
    """Current system status."""
    url = f'{config.get("JAWS", "url")}/status'
    try:
        r = requests.get(url)
    except requests.exceptions.RequestException:
        raise SystemExit("JAWS Central is DOWN")
    if r.status_code != 200:
        raise SystemExit(r.text)
    result = r.json()
    _print_json(result)


@main.command()
def info() -> None:
    """JAWS version and info."""
    url = f'{config.get("JAWS", "url")}/info'
    try:
        r = requests.get(url)
    except requests.exceptions.RequestException:
        raise SystemExit("JAWS Central is DOWN")
    if r.status_code != 200:
        raise SystemExit(r.text)
    result = r.json()
    _print_json(result)


@main.command()
@click.argument("src_dir")
@click.argument("dest_dir")
@click.option("--flatten", is_flag=True, default=False, help="Flatten shard dirs")
def wfcopy(src_dir: str, dest_dir: str, flatten) -> None:
    """Simplify Cromwell output."""
    wfc.wfcopy(src_dir, dest_dir, flatten)


def _request(rest_op, url, data={}, files={}):
    """Perform specified REST operation.  A JSON response is expected."""
    access_token = config.get("USER", "token")
    if not access_token:
        raise SystemExit("User access token required; contact an admin to get yours.")
    header = {"Authorization": f"Bearer {access_token}"}
    response = None
    try:
        if rest_op == "GET":
            response = requests.get(url, headers=header)
        elif rest_op == "PUT":
            response = requests.put(url, headers=header)
        elif rest_op == "POST":
            response = requests.post(url, data=data, files=files, headers=header)
        else:
            raise ValueError(f"Unsupported REST request type: {rest_op}")
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
    if response.status_code < 200 or response.status_code > 299:
        try:
            result = response.json()
        except Exception:
            raise SystemExit(response.text)
        if "error" in result:
            raise SystemExit(result["error"])
        else:
            raise SystemExit(result)
    return response.json()


def _print_json(j):
    click.echo(json.dumps(j, indent=4, sort_keys=True))


@main.command()
@click.option(
    "--site", default="ALL", help="limit results to this compute-site; default=all"
)
def queue(site: str) -> None:
    """List of user's current runs"""
    data = {
        "delta_days": 0,
        "site_id": site.upper(),
        "active_only": True,
        "result": "any",
    }
    url = f'{config.get("JAWS", "url")}/search'
    result = _request("POST", url, data)
    _print_json(result)


@main.command()
@click.option("--days", default=1, help="history going back this many days; default=1")
@click.option(
    "--site", default="ALL", help="limit results to this compute-site; default=all"
)
@click.option(
    "--result", default="any", help="limit results to this result; default=any"
)
def history(days: int, site: str, result: str) -> None:
    """Print a list of the user's past runs."""
    if days < 1:
        raise SystemExit("User error: --days must be a positive integer")
    if site:
        site = site.upper()
    data = {
        "delta_days": days,
        "site_id": site.upper(),
        "active_only": False,
        "result": result.lower(),
    }
    url = f'{config.get("JAWS", "url")}/search'
    result = _request("POST", url, data)
    _print_json(result)


def _run_status(run_id: int) -> Dict[str, str]:
    """Return the status of a run in JSON format."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}'
    return _request("GET", url)


@main.command()
@click.argument("run_id")
def status(run_id: int) -> None:
    """Print the current status of a run."""

    result = _run_status(run_id)
    _print_json(result)


@main.command()
@click.argument("run_id")
@click.option("--fmt", default="text", help="the desired output format: [text|json]")
def task_status(run_id: int, fmt: str) -> None:
    """Show the current status of each Task."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/task_status'
    result = _request("GET", url)
    if fmt == "json":
        _print_json(result)
    else:
        click.echo(
            "#CROMWELL_RUN_ID\tTASK_NAME\tATTEMPT\tCROMWELL_JOB_ID\tSTATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON"
        )
        for row in result:
            row[2] = str(row[2])
            row[3] = str(row[3])
            click.echo("\t".join(row))


@main.command()
@click.argument("run_id")
def metadata(run_id: int) -> None:
    """ Print detailed metadata for a run, generated by cromwell. """

    url = f'{config.get("JAWS", "url")}/run/{run_id}/metadata'
    result = _request("GET", url)
    _print_json(result)


@main.command()
@click.argument("run_id")
@click.option("--fmt", default="text", help="the desired output format: [text|json]")
def log(run_id: int, fmt: str) -> None:
    """View the log of Run state transitions for the workflow as a whole."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/run_log'
    result = _request("GET", url)
    if fmt == "json":
        _print_json(result)
    else:
        click.echo("#STATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON")
        for log_entry in result:
            click.echo("\t".join(log_entry))


@main.command()
@click.argument("run_id")
@click.option("--fmt", default="text", help="the desired output format: [text|json]")
def task_log(run_id: int, fmt: str) -> None:
    """Get log of each Task's state transitions."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/task_log'
    result = _request("GET", url)
    if fmt == "json":
        _print_json(result)
    else:
        click.echo(
            "#CROMWELL_RUN_ID\tTASK_NAME\tATTEMPT\tCROMWELL_JOB_ID\tSTATUS_FROM\tSTATUS_TO\tTIMESTAMP\tREASON"
        )
        for row in result:
            row[2] = str(row[2])
            row[3] = str(row[3])
            click.echo("\t".join(row))


@main.command()
@click.argument("run_id")
@click.option("--fmt", default="text", help="the desired output format: [text|json]")
def errors(run_id: int, fmt: str) -> None:
    """View error messages and stderr for failed Tasks."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/errors'
    result = _request("GET", url)
    if fmt == "json":
        _print_json(result)
    else:
        for task_name in result:
            click.echo(f"{task_name}:")
            click.echo(result[task_name])
            click.echo("\n")


@main.command()
@click.argument("run_id")
def cancel(run_id):
    """Cancel a run; prints whether aborting was successful or not."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/cancel'
    result = _request("PUT", url)
    _print_json(result)


@main.command()
def cancel_all():
    """Cancel all active runs."""
    url = f'{config.get("JAWS", "url")}/run/cancel-all'
    result = _request("PUT", url)
    _print_json(result)


def _list_sites() -> None:
    """List available Sites."""

    url = f'{config.get("JAWS", "url")}/site'
    result = _request("GET", url)
    _print_json(result)


@main.command()
def list_sites() -> None:
    """List available compute Sites"""
    _list_sites()


@main.command()
@click.argument("wdl_file", nargs=1)
@click.argument("json_file", nargs=1)
@click.argument("site", nargs=1)
@click.option("--tag", help="identifier for the run")
@click.option("--no-cache", is_flag=True, help="Disable call-caching for this run")
def submit(wdl_file: str, json_file: str, site: str, tag: str, no_cache: bool):
    """Submit a run for execution at a JAWS-Site.
    Available sites can be found by running 'jaws run list-sites'.
    """
    from jaws_client import workflow

    # the users' jaws id may not match the linux uid where the client is installed
    url = f'{config.get("JAWS", "url")}/user'
    result = _request("GET", url)
    uid = result["uid"]

    staging_subdir = config.get("JAWS", "staging_dir")
    staging_user_subdir = os.path.join(staging_subdir, uid)
    globus_host_path = config.get("GLOBUS", "host_path")
    output_directory = config.get("JAWS", "data_repo_basedir")
    input_site_id = config.get("JAWS", "site_id")
    local_staging_endpoint = workflow.join_path(globus_host_path, staging_user_subdir)

    # GET SITE INFO
    compute_site_id = site.upper()
    url = f'{config.get("JAWS", "url")}/site/{compute_site_id}'
    result = _request("GET", url)
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
        staged_wdl, zip_file = wdl.compress_wdls(local_staging_endpoint)
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

    # turning off call-caching requires a Cromwell options json file be created
    options_json_file = None
    if no_cache is True:
        options_json_file = workflow.join_path(
            local_staging_endpoint, f"{submission_id}.options.json"
        )
        with open(options_json_file, "w") as fh:
            fh.write('{"read_from_cache": false, "write_to_cache": false}')

    # write the file transfer manifest; jaws-central shall submit the transfer to globus
    manifest_file = workflow.Manifest(local_staging_endpoint, compute_uploads_subdir)
    manifest_file.add(staged_wdl, zip_file, staged_json, orig_json, *moved_files)
    if options_json_file:
        manifest_file.add(options_json_file)
    staged_manifest = workflow.join_path(staging_user_subdir, f"{submission_id}.tsv")
    manifest_file.write_to(staged_manifest)

    # SUBMIT RUN TO CENTRAL
    local_endpoint_id = config.get("GLOBUS", "endpoint_id")
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
    url = f'{config.get("JAWS", "url")}/run'
    logger.debug(f"Submitting run: {data}")
    result = _request("POST", url, data, files)
    result["max_ram_gb"] = max_ram_gb
    del result["output_dir"]
    _print_json(result)


@main.command()
@click.argument("wdl_file", nargs=1)
def inputs(wdl_file: str) -> None:
    """Generate inputs template (JSON) from workflow (WDL) file."""

    from jaws_client import workflow

    if not os.path.isfile(wdl_file):
        raise IOError(f"File not found: {wdl_file}")
    stdout, stderr = workflow.womtool("inputs", wdl_file)
    if stderr:
        raise SystemExit(stderr)
    click.echo(stdout.strip())


@main.command()
@click.argument("wdl_file", nargs=1)
def validate(wdl_file: str) -> None:
    """Validate a WDL using Cromwell's WOMTool."""

    from jaws_client import workflow

    if not os.path.isfile(wdl_file):
        raise IOError(f"File not found: {wdl_file}")
    stdout, stderr = workflow.womtool("inputs", wdl_file)
    if stderr:
        raise SystemExit(stderr)
    else:
        click.echo("Workflow is OK")


@main.command()
@click.argument("run_id")
@click.argument("dest")
def get(run_id: int, dest: str) -> None:
    """Copy the output of a run to the specified folder."""

    from jaws_client import workflow

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
            [
                "-rLtq",
                "--chmod=Du=rwx,Dg=rwx,Do=,Fu=rw,Fg=rw,Fo=",
                "--exclude=inputs",
                "--exclude=tmp.*",
            ],
        )
    except IOError as error:
        logger.error(f"Rsync output failed for run {run_id}: {error}")
        raise SystemExit(f"Error getting output for run {run_id}: {error}")
    if result.returncode != 0:
        err_msg = f"Failed to rsync {src}->{dest}: {result.stdout}; {result.stderr}"
        raise SystemExit(err_msg)


@main.command()
@click.argument("uid")
@click.argument("email")
@click.argument("name")
@click.option("--admin", is_flag=True, default=False, help="Grant admin access")
def add_user(uid: str, email: str, name: str, admin: bool) -> None:
    """Add new user and get JAWS OAuth access token (restricted)."""

    data = {"uid": uid, "email": email, "name": name, "admin": admin}
    url = f'{config.get("JAWS", "url")}/user'
    result = _request("POST", url, data)
    _print_json(result)


def jaws():
    """Entrypoint for jaws-client app."""
    main.add_command(deprecated.run)
    main()
