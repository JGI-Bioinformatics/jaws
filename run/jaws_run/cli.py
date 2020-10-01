#!/usr/bin/env python

"""
JAWS Run Service
"""

import os
import click

from jaws_run import log, config

JAWS_LOG_ENV = "JAWS_RUN_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), f"{__package__}.log")
JAWS_CONFIG_ENV = "JAWS_RUN_CONFIG"
JAWS_CWD_CONFIG = os.path.join(os.getcwd(), f"{__package__}.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(config_file: str, log_file: str, log_level: str):
    """JAWS Run Service"""
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
def rpc() -> None:
    """Start RPC server for Central."""
    from jaws_rpc import rpc_server
    from jaws_run.rpc_methods import rpc_methods

    rpc_server_params = config.conf.get_section("RPC")
    app = rpc_server.RpcServer(rpc_server_params, rpc_methods)
    app.start_server()


@cli.command()
def daemon() -> None:
    """Start daemon."""

    # create tables if not exists
    from jaws_run import db
    session = db.Session()
    db.create_all(db.engine, session)

    # start daemon process
    from jaws_run import daemon
    jawsd = daemon.Daemon()
    jawsd.start_daemon()


def jaws():
    """Entrypoint for jaws-run app."""
    cli()
