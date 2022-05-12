#!/usr/bin/env python3

"""JAWS Central CLI"""

import os
import click
import logging
from jaws_central import config, log
import connexion
from flask_cors import CORS
from urllib.parse import quote_plus
from sqlalchemy.pool import QueuePool
from jaws_central.models_fsa import db
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


@cli.command()
def auth() -> None:
    """Start JAWS OAuth server"""
    logger = logging.getLogger(__package__)
    logger.debug("Initializing OAuth server")
    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp(__name__, specification_dir=basedir)
    connex.add_api("swagger.auth.yml")

    connex.app.config["SQLALCHEMY_DATABASE_URI"] = "%s://%s:%s@%s:%s/%s" % (
        config.conf.get("DB", "dialect"),
        config.conf.get("DB", "user"),
        quote_plus(config.conf.get("DB", "password")),
        config.conf.get("DB", "host"),
        config.conf.get("DB", "port"),
        config.conf.get("DB", "db"),
    )
    connex.app.config["SQLALCHEMY_ECHO"] = False
    connex.app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False
    connex.app.config["SQLALCHEMY_ENGINE_OPTIONS"] = {
        "poolclass": QueuePool,
        "pool_pre_ping": True,
        "pool_recycle": 3600,
        "pool_size": 5,
        "max_overflow": 10,
        "pool_timeout": 30,
    }
    db.init_app(connex.app)

    # create tables if not exists
    with connex.app.app_context():
        try:
            db.create_all()
            db.session.commit()
        except Exception as error:
            db.session.rollback()
            logger.exception(f"Failed to create tables: {error}")
            raise

    # define port
    port = int(config.conf.get("HTTP", "auth_port"))  # defaults to 3000

    # start OAuth server
    connex.run(host="0.0.0.0", port=port, debug=False)


@cli.command()
def rest() -> None:
    """Start JAWS REST server."""
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

    connex.app.config["SQLALCHEMY_DATABASE_URI"] = "%s://%s:%s@%s:%s/%s" % (
        config.conf.get("DB", "dialect"),
        config.conf.get("DB", "user"),
        quote_plus(config.conf.get("DB", "password")),
        config.conf.get("DB", "host"),
        config.conf.get("DB", "port"),
        config.conf.get("DB", "db"),
    )
    connex.app.config["SQLALCHEMY_ECHO"] = False
    connex.app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False
    connex.app.config["SQLALCHEMY_ENGINE_OPTIONS"] = {
        "poolclass": QueuePool,
        "pool_pre_ping": True,
        "pool_recycle": 3600,
        "pool_size": 5,
        "max_overflow": 10,
        "pool_timeout": 30,
    }

    # add CORS support
    CORS(connex.app)

    db.init_app(connex.app)

    # create tables if not exists
    with connex.app.app_context():
        try:
            db.create_all()
            db.session.commit()
        except Exception as error:
            db.session.rollback()
            logger.exception(f"Failed to create tables: {error}")
            raise

    # init RPC clients
    site_rpc_params = config.conf.get_all_sites_rpc_params()
    rpc_index.rpc_index = rpc_index.RpcIndex(site_rpc_params, logger)

    # define port
    port = int(config.conf.get("HTTP", "rest_port"))  # defaults to 5000

    # start REST server
    connex.run(host="0.0.0.0", port=port, debug=False)


@cli.command()
def rpc() -> None:
    """Start JAWS-Central RPC server."""
    from jaws_rpc import rpc_server
    from jaws_central.database import Session
    from jaws_central import rpc_operations

    rpc_params = config.conf.get_section("RPC_SERVER")
    logger = logging.getLogger(__package__)
    app = rpc_server.RpcServer(rpc_params, logger, rpc_operations.operations, Session)
    app.start_server()


@cli.command()
def run_daemon() -> None:
    """Start run-daemon"""
    from jaws_central.database import engine, Session
    from jaws_central import models
    from jaws_central.run_daemon import RunDaemon

    session = Session()
    models.create_all(engine, session)

    logger = logging.getLogger(__package__)

    # init RPC clients
    site_rpc_params = config.conf.get_all_sites_rpc_params()
    rpc_index.rpc_index = rpc_index.RpcIndex(site_rpc_params, logger)

    #    import daemon
    #    import lockfile
    #    pidfile = os.path.join(config.conf.get("JAWS", "lock_dir"), "rund.pid")
    #    with daemon.DaemonContext(pidfile=lockfile.FileLock(pidfile)):
    #        rund = RunDaemon(rpc_index.rpc_index)
    #        rund.start_daemon()

    rund = RunDaemon(rpc_index.rpc_index)
    rund.start_daemon()


def jaws():
    """Entrypoint for jaws-server app."""
    cli()
