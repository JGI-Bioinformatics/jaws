#!/usr/bin/env python

"""
JAWS Task Service
"""

import os
import click
from jaws_task import log, config
from jaws_task.rpc_operations import operations
from jaws_rpc import rpc_server

JAWS_LOG_ENV = "JAWS_TASK_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), f"{__package__}.log")
JAWS_CONFIG_ENV = "JAWS_TASK_CONFIG"
JAWS_CWD_CONFIG = os.path.join(os.getcwd(), f"{__package__}.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(config_file: str, log_file: str, log_level: str):
    """JAWS-Site"""
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

    # create tables if not exists
    from jaws_task import db
    session = db.Session()
    db.create_all(db.engine, session)

@cli.command()
def server() -> None:
    """Start RPC server."""

    rpc_server_params = config.conf.get_section("RPC_SERVER")
    app = rpc_server.RpcServer(rpc_server_params, operations)
    app.start_server()


def jaws():
    """Entrypoint for jaws-site app."""
    cli()
