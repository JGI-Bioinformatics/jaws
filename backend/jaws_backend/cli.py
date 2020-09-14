#!/usr/bin/env python

"""
This script provides a CLI interface to the JAWS backend, to be used by Cromwell for submitting,
killing, and checking on the status of tasks via JSON-RPC vs jaws-site.
"""

import os
import click
import configparser
from jaws_rpc import rpc_client
from jaws_backend import log

JAWS_LOG_ENV = "JAWS_BACKEND_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), f"{__package__}.log")
JAWS_CONFIG_ENV = "JAWS_SITE_CONFIG"
JAWS_CWD_CONFIG = os.path.join(os.getcwd(), f"{__package__}.conf")

conf = configparser.ConfigParser()
logger = None


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(config_file: str, log_file: str, log_level: str):
    """JAWS Cromwell Backend"""
    global logger
    global conf
    if log_file is None:
        log_file = (
            os.environ[JAWS_LOG_ENV] if JAWS_LOG_ENV in os.environ else JAWS_CWD_LOG
        )
    logger = log.setup_logger(__package__, log_file, log_level)
    if config_file is None:
        config_file = (
            os.environ[JAWS_CONFIG_ENV]
            if JAWS_CONFIG_ENV in os.environ
            else JAWS_CWD_CONFIG
        )
    conf.read(config_file)
    if "SITE_RPC" not in conf:
        error = f"Config missing SITE_RPC section: {config_file}"
        logger.error(error)
        raise SystemExit(error)
    logger.info(f"Config using {config_file}")


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
def submit(script, job_name, cwd, out, err, max_time, memory_gb) -> None:
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
    global logger
    logger.debug(f"RPC {operation} with {params}")
    try:
        with rpc_client.RpcClient(conf["SITE_RPC"]) as rpc:
            response = rpc.request(operation, params)
    except Exception as error:
        logger.error(f"RPC {operation} failed: {error}")
        raise
    if "result" in response:
        logger.debug(f"- result: {response['result']}")
        print(response["result"])
    elif "error" in response:
        logger.error(
            f"RPC {operation} with {params} failed: {response['error']['message']}"
        )
        raise SystemExit(response["error"]["message"])
    else:
        logger.error(
            f"RPC {operation} with {params} returned invalid response: {response}"
        )
        raise SystemExit(f"Invalid JSON-RPC2 response: {response}")
