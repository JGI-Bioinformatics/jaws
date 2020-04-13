#!/usr/bin/env python

"""
JAWS Site server runs at each computing site and is comprised of:
(a) RPC server for handling user (sync) requests and
(b) daemon for performing periodic (async) maintenance tasks.
Each computing site also has a Cromwell server instance, typically installed on the same server.
"""

import os
import click

from jaws_site import database, rpc_server, jawsd, log, config

JAWS_SITE_LOG = "JAWS_SITE_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), 'jaws-site.log')


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
def cli(config_file: str, log_file: str):
    """JAWS-Site."""
    # Initialize logging and configuration singletons;
    # as they are singletons, the Click context object is not needed.
    if log_file is None:
        log_file = os.environ[JAWS_SITE_LOG] if JAWS_SITE_LOG in os.environ else JAWS_CWD_LOG
    logger = log.setup_logger(__package__, log_file)
    conf = config.Configuration(config_file)
    logger.debug(f"Config using {conf.config_file}")


@cli.command()
def server() -> None:
    """Start RPC server."""
    app = rpc_server.RpcServer()
    app.start_server()


@cli.command()
def daemon() -> None:
    """Start daemon."""
    db = database.Database(config.conf)
    daemon = jawsd.Daemon(db)
    daemon.start_daemon()


def jaws():
    """Entrypoint for jaws-site app."""
    cli()
