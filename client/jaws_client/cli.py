"""
JAWS CLI
"""

import sys
import click
import os
import pathlib
import requests
import json
import subprocess
import warnings
from typing import Dict
from jaws_client import log as logging
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
    except Exception as error:
        raise (f"Unable to communicate with JAWS Central: {error}")
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
    except Exception as error:
        raise (f"Unable to communicate with JAWS Central: {error}")
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
    except Exception as error:
        raise (f"Unable to communicate with JAWS Central: {error}")
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
    help="List runs from all users; default=False",
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
def outfiles(run_id: int) -> None:
    """Print list of a run's output files."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/outfiles'
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
    """Get Task logs and current status"""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/task_log'
    result = _request("GET", url)
    header = [
        "NAME",
        "CACHED",
        "STATUS",
        "QUEUED",
        "RUNNING",
        "FINISHED",
        "QUEUE_DUR",
        "RUN_DUR",
    ]
    _convert_all_fields_to_localtime(result, columns=[3, 4, 5])
    if fmt == "json":
        _print_json(result)
    elif fmt == "tab":
        _print_tab_delimited_table(header, result)
    else:
        _print_space_delimited_table(header, result)


@main.command()
@click.argument("run_id")
def task_summary(run_id: int) -> None:
    """Get summary info on each Task"""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/task_summary'
    result = _request("GET", url)
    _convert_all_fields_to_localtime(
        result, keys=["queue_start", "run_start", "run_end"]
    )
    _print_json(result)


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
def errors(run_id: int) -> None:
    """View error messages and stderr for failed Tasks."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/errors'
    errors_report = _request("GET", url)
    _print_json(errors_report)


@main.command()
@click.argument("run_id")
def running_tasks(run_id: int) -> None:
    """View information about running Tasks."""

    url = f'{config.get("JAWS", "url")}/run/{run_id}/running'
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


def _get_uid() -> str:
    """Get the current user's JAWS user ID"""
    url = f'{config.get("JAWS", "url")}/user'
    result = _request("GET", url)
    uid = result["uid"]
    return uid


@main.command()
@click.argument("wdl_file", nargs=1)
@click.argument("json_file", nargs=1)
@click.argument("site", nargs=1)
@click.option("--tag", help="identifier for the run")
@click.option("--no-cache", is_flag=True, help="Disable call-caching for this run")
@click.option("--quiet", is_flag=True, help="Don't print copy progress bar")
@click.option(
    "--sub", default=None, help="Subworkflows zip (optional; by default, auto-generate)"
)
@click.option(
    "--webhook",
    default=None,
    help="If provided, JAWS will POST to this URL when Run completes.",
)
def submit(
    wdl_file: str,
    json_file: str,
    site: str,
    tag: str,
    no_cache: bool,
    quiet: bool,
    sub: str,
    webhook: str,
):
    """Submit a run for execution at a JAWS-Site.
    Available sites can be found by running 'jaws run list-sites'.
    """
    from jaws_client import workflow
    from jaws_client.workflow import WdlError

    wdl_file = os.path.abspath(wdl_file)
    json_file = os.path.abspath(json_file)
    input_dir = config.get("JAWS", "inputs_dir")
    input_site_id = config.get("JAWS", "site_id")
    output_dir = config.get("JAWS", "downloads_dir")

    params = {}
    if sub:
        params["subworkflows_zip_file"] = sub
    if quiet:
        params["quiet"] = quiet
    try:
        run = workflow.Run(
            wdl_file, json_file, input_dir, input_site_id, output_dir, **params
        )
    except WdlError as error:
        raise SystemExit(f"Your workflow has an error: {error}")
    except Exception as error:
        raise SystemExit(f"Error initializing Run: {error}")

    # post Run to jaws-central
    data = {
        "input_site_id": input_site_id,
        "compute_site_id": site.upper(),
        "submission_id": run.submission_id,
        "caching": False if no_cache else True,
        "max_ram_gb": run.wdl.max_ram_gb,
        "manifest": json.dumps(run.manifest.files),
        # the following aren't used by JAWS but are recorded for users' benefit
        "wdl_file": wdl_file,
        "json_file": json_file,
        "tag": tag,
        "webhook": webhook,
    }
    url = f'{config.get("JAWS", "url")}/run'
    logger.debug(f"Submitting run: {data}")
    result = _request("POST", url, data)
    result["max_ram_gb"] = run.wdl.max_ram_gb
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

    run_info = _run_status(run_id, True)
    status = run_info["status"]
    submission_id = run_info["submission_id"]

    if status not in ["download complete", "email sent", "done"]:
        sys.exit(f"Run {run_id} output is not yet available; status is {status}")

    # this shouldn't be necessary as wdl, json, zip were copied to this folder during submission
    pathlib.Path(dest).mkdir(parents=True, exist_ok=True)

    src = None
    if run_info["input_site_id"] == run_info["compute_site_id"]:
        # get the results from the cromwell-execution dir
        url = f'{config.get("JAWS", "url")}/run/{run_id}/root'
        src = _request("GET", url)
    else:
        # results will be in the downloads dir
        src_base = config.get("JAWS", "downloads_dir")
        src = f"{src_base}/{submission_id}"

    # Write the "outputs.json" file because it contains non-file outputs (e.g. numbers)
    url = f'{config.get("JAWS", "url")}/run/{run_id}/outputs'
    outputs = _request("GET", url)
    if outputs is None:
        err_msg = f"There are no outputs for Run {run_id} at this time."
        sys.exit(err_msg)
    outputs_file = os.path.normpath(os.path.join(dest, "outputs.json"))
    if not quiet:
        print("Writing outputs.json")
    _write_json_file(outputs_file, outputs)

    # the original wdl, json, zip files were copied to the inputs dir during submission
    input_dir = config.get("JAWS", "inputs_dir")
    wdl_src_file = os.path.join(input_dir, f"{submission_id}.wdl")
    wdl_dest_file = os.path.join(dest, os.path.basename(run_info["wdl_file"]))
    _copy_file(wdl_src_file, wdl_dest_file, quiet)
    json_src_file = os.path.join(input_dir, f"{submission_id}.json")
    json_dest_file = os.path.join(dest, os.path.basename(run_info["json_file"]))
    _copy_file(json_src_file, json_dest_file, quiet)
    zip_src_file = os.path.join(input_dir, f"{submission_id}.zip")
    zip_dest_file = os.path.join(dest, "subworkflows.zip")
    if os.path.isfile(zip_src_file):
        _copy_file(zip_src_file, zip_dest_file, quiet)

    if complete is True:
        _get_complete(run_id, src, dest)
    else:
        _get_outfiles(run_id, src, dest, quiet)


def _write_json_file(outfile: str, contents: dict):
    try:
        with open(outfile, "w") as fh:
            fh.write(json.dumps(contents, sort_keys=True, indent=4))
        os.chmod(outfile, 0o0664)
    except Exception as error:
        sys.exit(f"Unable to write json file, {outfile}: {error}")


def rsync(src, dest, options=["-rLtq"]):
    """Copy source to destination using rsync.

    :param src: Source path
    :type src: str
    :param dest: Destination path
    :type dest: str
    :param options: rsync options
    :type options: str
    :return: None
    """
    return subprocess.run(
        ["rsync", *options, src, dest], capture_output=True, text=True
    )


def _get_complete(run_id: int, src: str, dest: str) -> None:
    """Copy the complete cromwell output dir"""
    src = f"{src}/"  # so rsync won't make an extra dir
    try:
        result = rsync(
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


def _get_outfiles(run_id: int, src_dir: str, dest_dir: str, quiet: bool) -> None:
    """Copy workflow outputs"""
    url = f'{config.get("JAWS", "url")}/run/{run_id}/outfiles'
    outfiles = _request("GET", url)
    if not outfiles or len(outfiles) == 0:
        err_msg = f"There are no outputs for Run {run_id} at this time."
        sys.exit(err_msg)
    for rel_path in outfiles:
        src_file = os.path.normpath(os.path.join(src_dir, rel_path))
        dest_file = os.path.normpath(os.path.join(dest_dir, rel_path))
        _copy_file(src_file, dest_file, quiet)


def _copy_file(src_file, dest_file, quiet=False):
    """Copy one output if it is a file"""
    try:
        dest_dir = os.path.dirname(dest_file)
        os.makedirs(dest_dir, exist_ok=True)
        copy_with_progress_bar(src_file, dest_file, quiet=quiet)
        os.chmod(dest_file, 0o0664)
    except IOError as error:
        warnings.warn(f"Error copying file, {src_file}: {error}")


def _convert_all_table_fields_to_localtime(table, **kwargs):
    for row in table:
        _convert_all_fields_to_localtime(row, **kwargs)


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
    local_datetime_obj = datetime_obj.replace(tzinfo=timezone.utc).astimezone(
        tz=local_tz_obj
    )
    return local_datetime_obj.strftime(fmt)


@main.command()
@click.argument("uid")
@click.argument("email")
@click.argument("name")
@click.argument("group")
@click.option("--admin", is_flag=True, default=False, help="Grant admin access")
def add_user(uid: str, email: str, name: str, group: str, admin: bool) -> None:
    """Add new user and get JAWS OAuth access token (restricted)."""

    data = {
        "uid": uid,
        "email": email,
        "name": name,
        "user_group": group,
        "admin": admin,
    }
    url = f'{config.get("JAWS", "url")}/user'
    result = _request("POST", url, data)
    _print_json(result)


def jaws():
    """Entrypoint for jaws-client app."""
    main()
