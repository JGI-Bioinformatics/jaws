#!/usr/bin/env python

"""
JAWS Site server runs at each computing site and is comprised of:
(a) RPC server for handling user (sync) requests and
(b) daemon for performing periodic (async) maintenance tasks.
Each computing site also has a Cromwell server instance, typically installed on the same server.
"""

import os
import click

from jaws_site import log, config

JAWS_LOG_ENV = "JAWS_SITE_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), f"{__package__}.log")
JAWS_CONFIG_ENV = "JAWS_SITE_CONFIG"
JAWS_CWD_CONFIG = os.path.join(os.getcwd(), f"{__package__}.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
def cli(config_file: str, log_file: str):
    """JAWS-Site."""
    # Initialize logging and configuration singletons;
    # as they are singletons, the Click context object is not needed.
    if log_file is None:
        log_file = os.environ[JAWS_LOG_ENV] if JAWS_LOG_ENV in os.environ else JAWS_CWD_LOG
    logger = log.setup_logger(__package__, log_file)
    if config_file is None:
        config_file = os.environ[JAWS_CONFIG_ENV] if JAWS_CONFIG_ENV in os.environ else JAWS_CWD_CONFIG
    conf = config.Configuration(config_file)
    if conf:
        logger.info(f"Config using {config_file}")


@cli.command()
def rpc() -> None:
    """Start RPC server."""
    from jaws_site import analysis
    from jaws_rpc import rpc_server

    site_rpc_server_params = config.conf.get_section("SITE_RPC_SERVER")
    app = rpc_server.RpcServer(site_rpc_server_params, analysis.operations)
    app.start_server()


@cli.command()
def daemon() -> None:
    """Start daemon."""
    from jaws_site import daemon

    jawsd = daemon.Daemon()
    jawsd.start_daemon()


def jaws():
    """Entrypoint for jaws-site app."""
    cli()
