#!/usr/bin/env python

"""
JAWS Site server runs at each computing site and is comprised of:
(a) RPC server for handling user (sync) requests and
(b) daemon for performing periodic (async) maintenance tasks.
Each computing site also has a Cromwell server instance, typically installed on the same server.
"""

import logging
import os
import sys
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
@click.option("--env-override",
              envvar="ENV_OVERRIDE_PREFIX",
              default=None,
              help="prefix for environment variable override of configuration values")
def cli(config_file: str, log_file: str, log_level: str, env_override: str):
    """JAWS-Site"""
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
    conf = config.Configuration(config_file, env_override)
    if conf:
        logger.info(f"Config using {config_file}")
    else:
        sys.exit("Unable to find config file.")

    from jaws_site.database import engine
    from jaws_site import models

    try:
        models.Base.metadata.create_all(bind=engine)
    except Exception as error:
        logger.exception(f"Failed to create db tables: {error}")


@cli.command()
def rpc_server() -> None:
    """Start RPC server."""
    # database must be imported after config
    from jaws_site.database import session_factory
    from jaws_site import rpc_operations
    from jaws_rpc import rpc_server
    from sqlalchemy.orm import scoped_session

    # start RPC server
    rpc_server_params = config.conf.get_section("RMQ")
    rpc_server_params["queue"] = config.conf.get("SITE", "id")
    logger = logging.getLogger(__package__)
    app = rpc_server.RpcServer(
        rpc_server_params, logger, rpc_operations.operations, scoped_session(session_factory)
    )
    app.start_server()


@cli.command()
def run_daemon() -> None:
    """Start run daemon."""
    from jaws_site.run_daemon import RunDaemon

    rund = RunDaemon()
    rund.start_daemon()


@cli.command()
def transfer_daemon() -> None:
    """Start transfer daemon."""
    # database must be imported after config
    from jaws_site.transfer_daemon import TransferDaemon

    transferd = TransferDaemon()
    transferd.start_daemon()


@cli.command()
def perf_metrics_daemon() -> None:
    """Start performance metrics daemon"""
    from jaws_site.perf_metrics_daemon import PerformanceMetricsDaemon

    perfd = PerformanceMetricsDaemon()
    perfd.start_daemon()


@cli.command()
def pool_manager_daemon() -> None:
    """Start pool_manager daemon."""
    from jaws_site.pool_manager_daemon import PoolManagerDaemon

    pool_managerd = PoolManagerDaemon()
    pool_managerd.start_daemon()


def jaws():
    """Entrypoint for jaws-site app."""
    cli()
