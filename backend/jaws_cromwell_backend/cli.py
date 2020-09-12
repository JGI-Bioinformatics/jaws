#!/usr/bin/env python

"""
This script provides a CLI interface to the JAWS backend, to be used by Cromwell for submitting,
killing, and checking on the status of tasks via JSON-RPC.
"""

import click
import configparser
from jaws_rpc import rpc_client

JAWS_USER_LOG = os.path.expanduser("~/jaws-backend.log")
JAWS_CONFIG_ENV = "JAWS_BACKEND_CONFIG"
JAWS_CONFIG_DEFAULT_FILE = os.path.expanduser("~/jaws-backend.conf")

conf = configparser.ConfigParser()


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config file")
def cli(config_file: str, log_file: str, log_level: str):
    """JAWS Cromwell Backend"""
    global conf
    if config_file is None:
        config_file = (
            os.environ[JAWS_CONFIG_ENV]
            if JAWS_CONFIG_ENV in os.environ
            else JAWS_CONFIG_DEFAULT_FILE
        )
    conf = config.Configuration(config_file)
    if "SITE_RPC" not in conf:
        raise SystemExit(f"Config file, {config_file}, is missing 'SITE_RPC'")


@cli.command()
@click.option("--script", "script", required=True, help="Script to run")
@click.option("--job-name", "job_name", required=True, help="Label to apply to job")
@click.option("--cwd", "cwd", required=True, help="Current working directory")
@click.option("--out", "out", required=True, help="File for standard output")
@click.option("--err", "err", required=True, help="File for standard error")
@click.option(
    "--max-time", "max_time", required=True, help="Maximum time as 'HH:MM:SS'"
)
@click.option(
    "--memory-gb", "memory_gb", required=True, help="Maximum RAM in gigabytes"
)
def submit() -> None:
    """Submit a job"""
    params = {
        "script": script,
        "job_name": job_name,
        "cwd": cwd,
        "out": out,
        "err": err,
        "max_time": max_time,
        "memory_gb": memory_gb,
    }
    __rpc("submit", params)


@cli.command()
@click.argument("job_id")
def kill(job_id: int) -> None:
    __rpc("kill", {"job_id": job_id})


@cli.command()
@click.argument("job_id")
def check_alive(job_id: int, flatten) -> None:
    __rpc("check_alive", {"job_id": job_id})


def jaws():
    """Entrypoint for jaws-cromwell-backend app."""
    cli()


def __rpc(operation: str, params: dict) -> None:
    """
    Perform RPC call.  If successful, print result, else raise exception.
    """
    global conf
    try:
        with rpc_client.RPC_Client(conf["SITE_RPC"]) as rpc:
            response = rpc.request(operation, params)
    except Exception as error:
        logger.error(f"RPC {operation} failed: {error}")
        raise
    if "result" in response:
        print(response["result"])
    elif "error" in response:
        raise SystemExit(response["error"]["message"])
    else:
        raise SystemExit(f"Invalid JSON-RPC2 response: {response}")
