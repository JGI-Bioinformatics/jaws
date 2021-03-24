"""
JAWS CLI
"""

import click
import os
import requests
import json
from jaws_client import log, config, analysis, catalog, admin
from jaws_client import wfcopy as wfc

JAWS_LOG_ENV = "JAWS_CLIENT_LOG"
JAWS_USER_LOG = os.path.expanduser("~/jaws.log")
JAWS_CONFIG_ENV = "JAWS_CLIENT_CONFIG"
JAWS_CONFIG_DEFAULT_FILE = os.path.expanduser("~/jaws-client.conf")
JAWS_USER_CONFIG_ENV = "JAWS_USER_CONFIG"
JAWS_USER_CONFIG_DEFAULT_FILE = os.path.expanduser("~/jaws.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "jaws_config_file", default=None, help="JAWS config file")
@click.option("--user", "user_config_file", default=None, help="User config file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(jaws_config_file: str, user_config_file: str, log_file: str, log_level: str):
    """JGI Analysis Workflows Service"""
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
    if user_config_file is None:
        user_config_file = (
            os.environ[JAWS_USER_CONFIG_ENV]
            if JAWS_USER_CONFIG_ENV in os.environ
            else JAWS_USER_CONFIG_DEFAULT_FILE
        )
    conf = config.Configuration(jaws_config_file, user_config_file)
    if conf:
        logger.debug(f"Config using {jaws_config_file}, {user_config_file}")


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
def info() -> None:
    """JAWS version and info."""
    url = f'{config.conf.get("JAWS", "url")}/info'
    try:
        r = requests.get(url)
    except requests.exceptions.RequestException:
        raise SystemExit("JAWS Central is DOWN")
    if r.status_code != 200:
        raise SystemExit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@cli.command()
@click.argument("src_dir")
@click.argument("dest_dir")
@click.option("--flatten", is_flag=True, default=False, help="Flatten shard dirs")
def wfcopy(src_dir: str, dest_dir: str, flatten) -> None:
    """Simplify Cromwell output."""
    wfc.wfcopy(src_dir, dest_dir, flatten)


def jaws():
    """Entrypoint for jaws-client app."""
    cli.add_command(analysis.run)
    cli.add_command(catalog.wdl)
    cli.add_command(admin.admin)
    cli()
