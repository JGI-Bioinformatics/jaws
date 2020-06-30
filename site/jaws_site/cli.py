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
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(config_file: str, log_file: str, log_level: str):
    """JAWS-Site.

    :param config_file: filename of configuration file
    :type config_file: str
    :param log_file: filename of log file
    :type log_file: str
    :param log_level: logging level
    :type log_level: str
    :return:
    """
    # Initialize logging and configuration singletons;
    # as they are singletons, the Click context object is not needed.
    if log_file is None:
        log_file = os.environ[JAWS_LOG_ENV] if JAWS_LOG_ENV in os.environ else JAWS_CWD_LOG
    logger = log.setup_logger(__package__, log_file, log_level)
    if config_file is None:
        config_file = os.environ[JAWS_CONFIG_ENV] if JAWS_CONFIG_ENV in os.environ else JAWS_CWD_CONFIG
    conf = config.Configuration(config_file)
    if conf:
        logger.info(f"Config using {config_file}")


@cli.command()
def central_rpc() -> None:
    """Start RPC server for Central."""
    from jaws_site import analysis
    from jaws_rpc import rpc_server

    central_rpc_server_params = config.conf.get_section("CENTRAL_RPC_SERVER")
    app = rpc_server.RpcServer(central_rpc_server_params, analysis.operations)
    app.start_server()


@cli.command()
def jtm_rpc() -> None:
    """Start RPC server for JTM."""
    from jaws_site import rpc_operations
    from jaws_rpc import rpc_server

    jtm_rpc_server_params = config.conf.get_section("LOCAL_RPC_SERVER")
    app = rpc_server.RpcServer(jtm_rpc_server_params, rpc_operations.operations)
    app.start_server()


@cli.command()
def daemon() -> None:
    """Start daemon."""

    # create tables if not exists
    from jaws_site.database import engine, Session
    from jaws_site import models
    session = Session()
    models.create_all(engine, session)

    # start daemon process
    from jaws_site import daemon
    jawsd = daemon.Daemon()
    jawsd.start_daemon()


def jaws():
    """Entrypoint for jaws-site app."""
    cli()
