#!/usr/bin/env python

"""
JAWS Site server runs at each computing site and is comprised of:
(a) RPC server for handling user (sync) requests and
(b) daemon for performing periodic (async) maintenance tasks.
Each computing site also has a Cromwell server instance, typically installed on the same server.
"""

import os
import sys
import click

from jaws_condor import log, config

JAWS_LOG_ENV = "JAWS_CONDOR_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), f"{__package__}.log")
JAWS_CONFIG_ENV = "JAWS_CONDOR_CONFIG"
JAWS_CWD_CONFIG = os.path.join(os.getcwd(), f"{__package__}.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(config_file: str, log_file: str, log_level: str):
    """JAWS-HTCondor"""
    # Initialize logging and configuration singletons;
    # as they are singletons, the Click context object is not needed.
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
    conf = config.Configuration(config_file)
    if conf:
        logger.info(f"Config using {config_file}")
    else:
        sys.exit("Unable to find config file.")


@cli.command()
def pool_manager_daemon() -> None:
    """Start pool_manager daemon."""

    from jaws_condor.pool_manager_daemon import PoolManagerDaemon

    pool_managerd = PoolManagerDaemon()
    pool_managerd.start_daemon()


def jaws():
    """Entrypoint for jaws-condor app."""
    cli()
