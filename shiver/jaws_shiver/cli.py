#!/usr/bin/env python

import click
import configparser
from jaws_rpc import rpc_client

JAWS_LOG_ENV = "JAWS_WORKFORCE_LOG"
JAWS_USER_LOG = os.path.expanduser("~/jaws-workforce.log")
JAWS_CONFIG_ENV = "JAWS_WORKFORCE_CONFIG"
JAWS_CONFIG_DEFAULT_FILE = os.path.expanduser("~/jaws-workforce.conf")

@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(config_file: str, log_file: str, log_level: str):
    """JAWS Workforce"""
    if log_file is None:
        log_file = (
            os.environ[JAWS_LOG_ENV] if JAWS_LOG_ENV in os.environ else JAWS_USER_LOG
        )
    logger = log.setup_logger(__package__, log_file, log_level)
    if jaws_config_file is None:
        jaws_config_file = (
            os.environ[JAWS_CONFIG_ENV]
            if JAWS_CONFIG_ENV in os.environ
            else JAWS_CONFIG_DEFAULT_FILE
        )
    conf = config.Configuration(config_file)
    if conf:
        logger.debug(f"Config using {config_file}")


@cli.command()
def status() -> None:
    """Current system status."""
    url = f'{config.conf.get("JAWS", "url")}/status'
    try:
        r = requests.get(url)
    except requests.exceptions.RequestException:
        raise SystemExit("JAWS Central is DOWN")
    if r.status_code != 200:
        raise SystemExit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@cli.command()
@click.argument("script")
def submit() -> None:
    """Submit a job"""
    response = _rpc('submit', {"script":script})


@cli.command()
@click.argument("job_id")
def cancel(job_id: int, flatten) -> None:
    response = _rpc('cancel', {"job_id":job_id})

@cli.command()
@click.argument("job_id")
def status(job_id: int, flatten) -> None:
    response = _rpc('status', {"job_id":job_id})

def jaws():
    """Entrypoint for jaws-workforce app."""
    #cli.add_command(accounting.run)
    cli()


def _rpc(operation, params):
    try:
        with rpc_client.RPC_Client(rpc_client_config) as rpc:
            response = rpc.request(operation, params)
            if "error" in response
                raise SystemExit(response["error"]["message"])
    except Exception as error:
        logger.error(f"RPC {operation} failed: {error}")
        raise
    if "result" in response:
        logger.info(f"Success: {response['result']}")
    else:
        logger.error(f"Failure: {response['error']['message']}")
        raise
