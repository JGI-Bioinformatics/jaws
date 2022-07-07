#!/usr/bin/env python3

"""JAWS Central CLI"""

import os
import click
import logging
from jaws_central import config, log
import connexion
from flask import _app_ctx_stack
from flask_cors import CORS
from sqlalchemy.orm import scoped_session
from jaws_rpc import rpc_index


JAWS_LOG_ENV = "JAWS_CENTRAL_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), f"{__package__}.log")
JAWS_CONFIG_ENV = "JAWS_CENTRAL_CONFIG"
JAWS_CWD_CONFIG = os.path.join(os.getcwd(), f"{__package__}.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(config_file: str, log_file: str, log_level: str) -> None:
    """JAWS-Central"""
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
        logger.debug(f"Config using {config_file}")

    from jaws_central.database import engine
    from jaws_central import models
    try:
        models.Base.metadata.create_all(bind=engine)
    except Exception as error:
        logger.exception(f"Failed to create db tables: {error}")


@cli.command()
def auth() -> None:
    """Start JAWS OAuth server"""
    # database must be initialized after config
    from jaws_central.database import session_factory

    logger = logging.getLogger(__package__)
    logger.debug("Initializing OAuth server")
    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp(__name__, specification_dir=basedir)
    connex.add_api("swagger.auth.yml")
    CORS(connex.app)
    connex.app.session = scoped_session(session_factory, scopefunc=_app_ctx_stack)

    @connex.app.teardown_appcontext
    def remove_session(*args, **kwargs):
        connex.app.session.remove()

    port = int(config.conf.get("HTTP", "auth_port"))  # defaults to 3000
    connex.run(host="0.0.0.0", port=port, debug=False)


@cli.command()
def rest() -> None:
    """Start JAWS REST server."""
    # database must be initialized after config
    from jaws_central.database import session_factory

    logger = logging.getLogger(__package__)
    logger.debug("Starting jaws-central REST server")
    auth_url = config.conf.get("HTTP", "auth_url")
    auth_port = config.conf.get("HTTP", "auth_port")
    if not auth_url.startswith("http"):
        auth_url = f"http://{auth_url}"
    os.environ["TOKENINFO_URL"] = f"{auth_url}:{auth_port}/tokeninfo"
    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp("JAWS_REST", specification_dir=basedir)
    connex.add_api("swagger.rest.yml")
    CORS(connex.app)
    connex.app.session = scoped_session(session_factory, scopefunc=_app_ctx_stack)

    @connex.app.teardown_appcontext
    def remove_session(*args, **kwargs):
        connex.app.session.remove()

    # init RPC clients
    site_rpc_params = config.conf.get_all_sites_rpc_params()
    rpc_index.rpc_index = rpc_index.RpcIndex(site_rpc_params, logger)

    port = int(config.conf.get("HTTP", "rest_port"))  # defaults to 5000
    connex.run(host="0.0.0.0", port=port, debug=False)


@cli.command()
def rpc() -> None:
    """Start JAWS-Central RPC server."""
    # database must be initialized after config
    from jaws_central.database import session_factory
    from jaws_rpc import rpc_server
    from jaws_central import rpc_operations

    rpc_params = config.conf.get_section("RPC_SERVER")
    logger = logging.getLogger(__package__)
    app = rpc_server.RpcServer(rpc_params, logger, rpc_operations.operations, scoped_session(session_factory))
    app.start_server()


@cli.command()
def run_daemon() -> None:
    """Start run-daemon"""
    from jaws_central.run_daemon import RunDaemon

    logger = logging.getLogger(__package__)

    # init RPC clients
    site_rpc_params = config.conf.get_all_sites_rpc_params()
    rpc_index.rpc_index = rpc_index.RpcIndex(site_rpc_params, logger)

    rund = RunDaemon(rpc_index.rpc_index)
    rund.start_daemon()


def jaws():
    """Entrypoint for jaws-server app."""
    cli()
