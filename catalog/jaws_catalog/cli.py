#!/usr/bin/env python3

"""JAWS Catalog CLI"""

import os
import click
import logging
import connexion
from urllib.parse import quote_plus
from jaws_catalog import config, log
from jaws_catalog.db import db


JAWS_LOG_ENV = "JAWS_CATALOG_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), f"{__package__}.log")
JAWS_CONFIG_ENV = "JAWS_CATALOG_CONFIG"
JAWS_CWD_CONFIG = os.path.join(os.getcwd(), f"{__package__}.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(config_file: str, log_file: str, log_level: str) -> None:
    """JAWS-Catalog"""
    # Initialize logging and configuration singletons;
    # as they are singletons, the Click context object is not needed.
    if log_file is None:
        log_file = os.environ[JAWS_LOG_ENV] if JAWS_LOG_ENV in os.environ else JAWS_CWD_LOG
    logger = log.setup_logger(__package__, log_file, log_level)
    if config_file is None:
        config_file = os.environ[JAWS_CONFIG_ENV] if JAWS_CONFIG_ENV in os.environ else JAWS_CWD_CONFIG
    conf = config.Configuration(config_file)
    if conf:
        logger.debug(f"Config using {config_file}")


@cli.command()
def server() -> None:
    """Start JAWS Catalog REST server."""
    logger = logging.getLogger(__package__)
    logger.debug("Starting jaws-catalog REST server")
    auth_url = config.conf.get("HTTP", "auth_url")
    auth_port = config.conf.get("HTTP", "auth_port")
    if not auth_url.startswith("http"):
        auth_url = f"http://{auth_url}"
    os.environ["TOKENINFO_URL"] = f"{auth_url}:{auth_port}/tokeninfo"
    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp("JAWS_REST", specification_dir=basedir)
    connex.add_api("swagger.yml")

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
        "pool_pre_ping": True,
        "pool_recycle": 300,
    }
    db.init_app(connex.app)

    # create tables if not exists
    with connex.app.app_context():
        try:
            db.create_all()
            db.session.commit()
        except Exception as e:
            logger.exception(f"Failed to create tables: {e}")
            raise

    # define port
    port = int(config.conf.get("HTTP", "rest_port"))  # defaults to 5000

    # start REST server
    connex.run(host="0.0.0.0", port=port, debug=False)


def jaws():
    """Entrypoint for jaws-server app."""
    cli()
