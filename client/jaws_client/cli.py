"""
JAWS CLI
"""

import sys
import click
import os
import requests
import json
from jaws_client import log, config, analysis, catalog


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """JGI Analysis Workflows Service"""
    pass


def find_config_file(config_file, env_var):
    """Find config file by checking env var, cwd, and home, in that order.

    :param config_file: filename of configuration file
    :type config_file: str
    :param env_var: name of config file environment variable
    :type env_var: str
    :return: absolute path if found, None otherwise.
    :rtype: str
    """
    if env_var is os.environ:
        cf = os.environ[env_var]
        if os.path.isfile(cf):
            return cf

    cf = os.path.abspath(config_file)
    if os.path.isfile(cf):
        return cf

    cf = os.path.join(os.environ["HOME"], config_file)
    if os.path.isfile(cf):
        return cf
    return None


@cli.command()
def status():
    """Current system status."""
    url = f'{config.conf.get("JAWS", "url")}/status'
    try:
        r = requests.get(url)
    except Exception:
        sys.exit("JAWS Central is DOWN")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


def jaws():
    """Entrypoint for jaws-client app."""
    log_file = os.path.join(os.environ["HOME"], "jaws.log")
    logger = log.setup_logger(__package__, log_file)

    config_file = find_config_file("jaws_client.ini", "JAWS_CLIENT_CONFIG")
    if not config_file:
        logger.critical("Unable to find jaws config file")
    config.conf = config.JawsConfig(config_file)

    # cli = click.CommandCollection(sources=[analysis.run, catalog.wdl)
    cli.add_command(analysis.run)
    cli.add_command(catalog.wdl)
    cli()
