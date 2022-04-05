"""
JAWS CLI
"""

import sys
import click
import os
import requests
import json
import uuid
import shutil
from typing import Dict
from jaws_client import log as logging
from jaws_client import deprecated
from jaws_client.config import Configuration
from jaws_client.copy_progress import copy_with_progress_bar

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
        sys.exit("JAWS Central is DOWN")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    _print_json(result)


@main.command()
def info() -> None:
    """JAWS version and info."""
    url = f'{config.get("JAWS", "url")}/info'
    try:
        r = requests.get(url)
    except requests.exceptions.RequestException:
        sys.exit("JAWS Central is DOWN")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    _print_json(result)


def _request(rest_op, url, data={}, files={}):
    """Perform specified REST operation.  A JSON response is expected."""
    access_token = config.get("USER", "token")
    if not access_token:
        sys.exit("User access token required; contact an admin to get yours.")
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
            sys.exit(f"Unsupported REST request type: {rest_op}")
    except requests.exceptions.Timeout as err:
        sys.exit("Unable to communicate with JAWS server (timeout)", err)
    except requests.exceptions.TooManyRedirects as err:
        sys.exit(
            "Unable to communicate with JAWS server (too many redirects; bad url?)", err
        )
    except requests.exceptions.HTTPError as err:
        sys.exit("Unable to communicate with JAWS server (http error)", err)
    except requests.exceptions.RequestException as err:
        sys.exit("Unable to communicate with JAWS server", err)
    if response.status_code < 200 or response.status_code > 299:
        try:
            result = response.json()
        except Exception:
            sys.exit(response.text)
        if "error" in result:
            sys.exit(result["error"])
        else:
            sys.exit(result)
    return response.json()


def _print_json(j):
    click.echo(json.dumps(j, indent=4, sort_keys=True))


@main.command()
@click.option(
    "--site", default="ALL", help="limit results to this compute-site; default=all"
)
@click.option(
    "--all", is_flag=True, default=False, help="List runs from all users; default=False"
)
def queue(site: str, all: bool) -> None:
    """List of user's current runs"""
    data = {
        "delta_days": 0,
        "site_id": site.upper(),
        "active_only": True,
        "result": "any",
        "all": all,
    }
    url = f'{config.get("JAWS", "url")}/search'
    result = _request("POST", url, data)
    _convert_all_fields_to_localtime(result, keys=["submitted", "updated"])
    _print_json(result)


@main.command()
@click.option("--days", default=1, help="history going back this many days; default=1")
@click.option(
    "--site", default="ALL", help="limit results to this compute-site; default=all"
)
@click.option(
    "--result", default="any", help="limit results to this result; default=any"
)
@click.option(
    "--all",
    is_flag=True,
    default=False,
    help="List runs from all users;d default=False",
)
def history(days: int, site: str, result: str, all: bool) -> None:
    """Print a list of the user's past runs."""
    if days < 1:
        sys.exit("User error: --days must be a positive integer")
    if site:
        site = site.upper()
    data = {
        "delta_days": days,
        "site_id": site.upper(),
        "active_only": False,
        "result": result.lower(),
        "all": all,
    }
    url = f'{config.get("JAWS", "url")}/search'
    result = _request("POST", url, data)
    _convert_all_fields_to_localtime(result, keys=["submitted", "updated"])
    _print_json(result)


def _run_status(run_id: int, verbose=False) -> Dict[str, str]:
    """Return the status of a run in JSON format."""

    if verbose is True:
        url = f'{config.get("JAWS", "url")}/run/{run_id}/complete'
    else:
        url = f'{config.get("JAWS", "url")}/run/{run_id}'
    return _request("GET", url)


@main.command()
@click.argument("run_id")
@click.option("--verbose", is_flag=True, help="Return all fields")
def status(run_id: int, verbose: bool) -> None:
    """Print the current status of a run."""

    result = _run_status(run_id, verbose)
    _convert_all_fields_to_localtime([result], keys=["submitted", "updated"])
    _print_json(result)


@main.command()
@click.argument("run_id")
@click.option("--fmt", default="text", help="the desired output format: [text|json]")
def task_status(run_id: int, fmt: str) -> None:
    """Show the current status of each Task."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/task_status'
    result = _request("GET", url)
    _convert_all_fields_to_localtime(result, columns=[4])
    header = [
        "NAME",
        "CROMWELL_JOB_ID",
        "CACHED",
        "STATUS",
        "TIMESTAMP",
        "COMMENT",
    ]
    if fmt == "json":
        _print_json(result)
    elif fmt == "tab":
        _print_tab_delimited_table(header, result)
    else:
        _print_space_delimited_table(header, result)


@main.command()
@click.argument("run_id")
def metadata(run_id: int) -> None:
    """Print detailed metadata for a run, generated by cromwell."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/metadata'
    result = _request("GET", url)
    _print_json(result)


@main.command()
@click.argument("run_id")
def outputs(run_id: int) -> None:
    """Print outputs for a run, generated by cromwell."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/outputs'
    result = _request("GET", url)
    _print_json(result)


@main.command()
@click.argument("run_id")
@click.option(
    "--fmt", default="text", help="the desired output format: [text|json|tab]"
)
def log(run_id: int, fmt: str) -> None:
    """View the log of Run state transitions for the workflow as a whole."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/run_log'
    result = _request("GET", url)
    _convert_all_fields_to_localtime(result, columns=[2])
    header = ["STATUS_FROM", "STATUS_TO", "TIMESTAMP", "COMMENT"]
    if fmt == "json":
        _print_json(result)
    elif fmt == "tab":
        _print_tab_delimited_table(header, result)
    else:
        _print_space_delimited_table(header, result)


def _make_table_printable(table):
    """
    Convert all of a table's data values to printable strings.
    """
    for row in table:
        for index in range(len(row)):
            if isinstance(row[index], bool):
                if row[index] is True:
                    row[index] = "True"
                elif row[index] is False:
                    row[index] = "False"
            if isinstance(row[index], (int, float)):
                row[index] = str(row[index])
            if row[index] is None:
                row[index] = ""
    return table


def _print_tab_delimited_table(header, table):
    """
    Print table with tab-separated columns, which is machine-parseable.
    """
    assert header
    table = _make_table_printable(table)
    header[0] = f"#{header[0]}"
    for row in table:
        click.echo("\t".join(header))
        for row in table:
            click.echo("\t".join(row))


def _print_space_delimited_table(header, table):
    """
    Print table with variable number of spaces to line-up columns for human-readability.
    """
    assert header
    table = _make_table_printable(table)
    header[0] = f"#{header[0]}"
    table.insert(0, header)
    col_widths = []

    # Get the max length of element in every col and add padding of 2 spaces
    for index in range(len(header)):
        col_widths.append(max([len(str(log_entry[index])) for log_entry in table]) + 2)
    for log_entry in table:
        click.echo(
            "".join(
                str(cell).ljust(col_widths[col_index])
                for col_index, cell in enumerate(log_entry)
            )
        )


@main.command()
@click.argument("run_id")
@click.option("--fmt", default="text", help="the desired output format: [text|json]")
def task_log(run_id: int, fmt: str) -> None:
    """Get log of each Task's state transitions."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/task_log'
    result = _request("GET", url)
    header = [
        "NAME",
        "CROMWELL_JOB_ID",
        "CACHED",
        "STATUS_FROM",
        "STATUS_TO",
        "TIMESTAMP",
        "COMMENT",
    ]
    _convert_all_fields_to_localtime(result, columns=[5])
    if fmt == "json":
        _print_json(result)
    elif fmt == "tab":
        _print_tab_delimited_table(header, result)
    else:
        _print_space_delimited_table(header, result)


def _convert_to_table(header: list, jdoc: dict):
    """Convert list of dictionaries to list of lists, given list of keys"""
    table = []
    for rec in jdoc:
        row = []
        for key in header:
            value = rec.get(key.lc(), None)
            row.append(value)
        table.append(row)
    return table


@main.command()
@click.argument("run_id")
@click.option("--fmt", default="text", help="the desired output format: [text|json]")
def task_summary(run_id: int, fmt: str) -> None:
    """Get summary of each Task's state durations."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/task_summary'
    result = _request("GET", url)
    header = [
        "NAME",
        "CROMWELL_JOB_ID",
        "CACHED",
        "RESULT",
        "QUEUED",
        "QUEUE_WAIT",
        "RUNTIME",
        "MAX_TIME",
    ]
    if fmt == "json":
        _print_json(result)
    elif fmt == "tab":
        _print_tab_delimited_table(header, result)
    else:
        _print_space_delimited_table(header, result)


@main.command()
@click.argument("run_id")
def errors(run_id: int) -> None:
    """View error messages and stderr for failed Tasks."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/errors'
    errors_report = _request("GET", url)
    _print_json(errors_report)


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
@click.option("--quiet", is_flag=True, help="Don't print copy progress bar")
@click.option("--default-container", default=None, help="The default Docker container to use for Tasks")
@click.option("--sub", default=None, help="Subworkflows zip (optional; by default, auto-generate)")
def submit(
    wdl_file: str,
    json_file: str,
    site: str,
    tag: str,
    no_cache: bool,
    quiet: bool,
    default_container: str,
    sub: str
):
    """Submit a run for execution at a JAWS-Site.
    Available sites can be found by running 'jaws run list-sites'.
    """
    user_zip_file = sub
    from jaws_client import workflow
    from jaws_client.workflow import WdlError

    wdl_file = os.path.abspath(wdl_file)
    json_file = os.path.abspath(json_file)

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
        wdl.validate(compute_max_ram_gb)
        max_ram_gb = wdl.max_ram_gb
    except WdlError as error:
        sys.exit(error)

    # any and all subworkflow WDL files must be supplied to Cromwell in a single ZIP archive
    try:
        staged_wdl, zip_file = wdl.prepare_wdls(local_staging_endpoint, user_zip_file)
    except IOError as error:
        sys.exit(f"Unable to copy WDLs to inputs dir: {error}")

    # VALIDATE INPUTS JSON
    try:
        inputs_json = workflow.WorkflowInputs(
            json_file, submission_id, wdl_loc=wdl_file
        )
    except json.JSONDecodeError as error:
        sys.exit(f"Your file, {json_file}, is not a valid JSON file: {error}")

    staged_json = workflow.join_path(local_staging_endpoint, f"{submission_id}.json")
    site_subdir = workflow.join_path(local_staging_endpoint, input_site_id)
    jaws_site_staging_dir = workflow.join_path(compute_basedir, compute_uploads_subdir)
    jaws_site_staging_site_subdir = workflow.join_path(
        jaws_site_staging_dir, input_site_id
    )

    # copy infiles in inputs json to site's inputs dir so they may be read by jaws user and
    # transferred to the compute site via globus
    try:
        moved_files = inputs_json.move_input_files(site_subdir, quiet)
    except Exception as error:
        sys.exit(f"Unable to copy input files: {error}")

    # the paths in the inputs json file are changed to their paths at the compute site
    modified_json = inputs_json.prepend_paths_to_json(jaws_site_staging_site_subdir)
    modified_json.write_to(staged_json)

    # the original inputs json file is kept with the run submission for record-keeping only
    orig_json = workflow.join_path(local_staging_endpoint, f"{submission_id}.orig.json")
    try:
        shutil.copy(json_file, orig_json)
    except IOError as error:
        sys.exit(f"Error copying JSON to {orig_json}: {error}")
    try:
        os.chmod(orig_json, 0o0664)
    except PermissionError as error:
        sys.exit(f"Unable to chmod {orig_json}: {error}")

    # specify Cromwell workflow options
    if default_container is None:
        default_container = config.get("JAWS", "default_container")
    workflow_options = {"default_runtime_attributes": {"docker": default_container}}
    if no_cache is True:
        workflow_options["read_from_cache"] = False
        workflow_options["write_to_cache"] = False

    # write workflow options json file
    options_json_file = workflow.join_path(
        local_staging_endpoint, f"{submission_id}.options.json"
    )
    with open(options_json_file, "w") as options_fh:
        json.dump(workflow_options, options_fh, indent=4)

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
        sys.exit(f"File not found: {wdl_file}")
    stdout, stderr = workflow.womtool("inputs", wdl_file)
    if stderr:
        sys.exit(stderr)
    click.echo(stdout.strip())


@main.command()
@click.argument("wdl_file", nargs=1)
def validate(wdl_file: str) -> None:
    """Validate a WDL using Cromwell's WOMTool."""

    from jaws_client import workflow

    if not os.path.isfile(wdl_file):
        sys.exit(f"File not found: {wdl_file}")
    stdout, stderr = workflow.womtool("inputs", wdl_file)
    if stderr:
        sys.exit(stderr)
    else:
        click.echo("Workflow is OK")


@main.command()
@click.argument("run_id")
@click.argument("dest")
@click.option(
    "--complete", is_flag=True, default=False, help="Get complete cromwell output"
)
@click.option(
    "--quiet", is_flag=True, default=False, help="Don't print copy progress bar"
)
def get(run_id: int, dest: str, complete: bool, quiet: bool) -> None:
    """Copy the output of a run to the specified folder."""

    result = _run_status(run_id, True)
    status = result["status"]
    src = result["output_dir"]

    if status != "download complete":
        sys.exit(f"Run {run_id} output is not yet available; status is {status}")

    if src is None:
        logger.error(f"Run {run_id} doesn't have an output_dir defined")
        sys.exit(f"Run {run_id} doesn't have an output_dir defined")

    if os.path.exists(dest) and os.path.isfile(dest):
        sys.exit(f"Error destination path is a file: {dest}")
    os.makedirs(dest, exist_ok=True)

    if complete is True:
        _get_complete(run_id, src, dest)
    else:
        _get_outputs(run_id, src, dest, quiet)


def _get_complete(run_id: int, src: str, dest: str) -> None:
    """Copy the complete cromwell output dir"""
    from jaws_client import workflow

    src = f"{src}/"  # so rsync won't make an extra dir
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
        sys.exit(f"Error getting output for run {run_id}: {error}")
    if result.returncode != 0:
        err_msg = f"Failed to rsync {src}->{dest}: {result.stdout}; {result.stderr}"
        sys.exit(err_msg)


def _get_outputs(run_id: int, src_dir: str, dest_dir: str, quiet: bool) -> None:
    """Copy workflow outputs"""
    url = f'{config.get("JAWS", "url")}/run/{run_id}/outputs'
    outputs = _request("GET", url)
    if not outputs:
        err_msg = f"There are no outputs for Run {run_id} at this time."
        sys.exit(err_msg)

    # Write the "outputs.json" file because it contains non-file outputs (e.g. numbers)
    outputs_file = os.path.normpath(os.path.join(dest_dir, "outputs.json"))
    try:
        with open(outputs_file, "w") as fh:
            fh.write(json.dumps(outputs, sort_keys=True, indent=4))
        os.chmod(outputs_file, 0o0664)
    except Exception as error:
        sys.exit(f"Unable to write output file: {error}")

    # the paths of workflow output files are listed in the outputs_file
    outputs = {}
    for (key, value) in outputs.items():
        if type(value) is list:
            for an_output in value:
                _copy_outfile(an_output, src_dir, dest_dir, quiet)
        else:
            _copy_outfile(value, src_dir, dest_dir, quiet)


def _copy_outfile(rel_path, src_dir, dest_dir, quiet=False):
    """Copy one output if it is a file"""
    if rel_path is None:
        return
    src_file = os.path.normpath(os.path.join(src_dir, rel_path))
    if not os.path.isfile(src_file):
        # not all of a workflow's "outputs" are files (may be string or number)
        return
    dest_file = os.path.normpath(os.path.join(dest_dir, rel_path))
    a_dest_dir = os.path.dirname(dest_file)
    os.makedirs(a_dest_dir, exist_ok=True)
    copy_with_progress_bar(src_file, dest_file, quiet=quiet)
    os.chmod(dest_file, 0o0664)


def _convert_all_fields_to_localtime(rec, **kwargs):
    if not rec:
        return
    for single_rec in rec:
        if "columns" in kwargs:
            for index in kwargs["columns"]:
                if single_rec[index]:
                    single_rec[index] = _utc_to_local(single_rec[index])
        elif "keys" in kwargs:
            for key in kwargs["keys"]:
                if key in single_rec and single_rec[key] is not None:
                    single_rec[key] = _utc_to_local(single_rec[key])


def _utc_to_local(utc_datetime):
    """Convert UTC time to the local time zone. This should handle daylight savings.
    Param:: utc_datetime: a string of date and time "2021-07-06 11:15:17".
    """
    from datetime import datetime, timezone
    import pytz

    # The timezone can be overwritten with a environmental variable.
    # JAWS_TZ should be set to a timezone in a similar format to 'US/Pacific'
    local_tz = os.environ.get("JAWS_TZ", None)
    local_tz_obj = ""
    if local_tz is None:
        local_tz_obj = datetime.now().astimezone().tzinfo
    else:
        local_tz_obj = pytz.timezone(local_tz)

    fmt = "%Y-%m-%d %H:%M:%S"
    datetime_obj = datetime.strptime(utc_datetime, fmt)
    local_datetime_obj = datetime_obj.replace(tzinfo=timezone.utc).astimezone(tz=local_tz_obj)
    return local_datetime_obj.strftime(fmt)


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
