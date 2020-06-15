#!/usr/bin/env python3

"""JAWS Central CLI"""

import os
import click
import logging
import connexion
from urllib.parse import quote_plus
import secrets
from jaws_central import config, log
from jaws_central.models import db
from jaws_rpc import rpc_index


JAWS_LOG_ENV = "JAWS_CENTRAL_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), f"{__package__}.log")
JAWS_CONFIG_ENV = "JAWS_CENTRAL_CONFIG"
JAWS_CWD_CONFIG = os.path.join(os.getcwd(), f"{__package__}.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
def cli(config_file: str, log_file: str) -> None:
    """JAWS-Central.

    :param config_file: filename of configuration file
    :type config_file: str
    :param log_file: filename of log file
    :type log_file: str
    :return:
    """
    # Initialize logging and configuration singletons;
    # as they are singletons, the Click context object is not needed.
    if log_file is None:
        log_file = os.environ[JAWS_LOG_ENV] if JAWS_LOG_ENV in os.environ else JAWS_CWD_LOG
    logger = log.setup_logger(__package__, log_file)
    if config_file is None:
        config_file = os.environ[JAWS_CONFIG_ENV] if JAWS_CONFIG_ENV in os.environ else JAWS_CWD_CONFIG
    conf = config.Configuration(config_file)
    if conf:
        logger.debug(f"Config using {config_file}")


@cli.command()
@click.option("--port", default=3000, help="Port (default=3000)")
def auth(port: int) -> None:
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
        "pool_pre_ping": True,
        "pool_recycle": 300,
    }
    db.init_app(connex.app)

    connex.run(host="0.0.0.0", port=port, debug=False)


@cli.command()
@click.option("--port", default=5000, help="Port (default=5000)")
def rest(port: int) -> None:
    """Start JAWS REST server."""
    logger = logging.getLogger(__package__)
    logger.debug("Starting jaws-central REST server")
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
        "pool_pre_ping": True,
        "pool_recycle": 300,
    }
    db.init_app(connex.app)

    rpc_index.rpc_index = rpc_index.RPC_Index(config.conf)

    connex.run(host="0.0.0.0", port=port, debug=False)


@cli.command()
def create_tables() -> None:
    """Create all database tables"""
    logger = logging.getLogger(__package__)
    logger.debug("Creating RDb tables")

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
        "pool_pre_ping": True,
        "pool_recycle": 300,
    }
    db.init_app(connex.app)

    with connex.app.app_context():
        logger.info("Creating all db tables")
        try:
            db.create_all()
            db.session.commit()
        except Exception as e:
            logger.exception(f"Failed to create tables: {e}")
            raise


@cli.command()
@click.argument("uid")
@click.argument("email")
@click.option("--admin", is_flag=True, default=False, help="Grant admin privileges")
def add_user(
    uid: str, email: str, admin: bool = False
) -> None:
    """Add user and generate OAuth2 token."""
    logger = logging.getLogger(__package__)
    logger.debug(f"Adding new user, {uid}")

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
        "pool_pre_ping": True,
        "pool_recycle": 300,
    }
    db.init_app(connex.app)

    with connex.app.app_context():
        # CHECK IF UID EXISTS
        from jaws_central.models import User
        user = db.session.query(User).get(uid)
        if user is not None:
            msg = f"Cannot add user {uid}; user.id already taken."
            logger.debug(msg)
            raise ValueError(msg)

        # GENERATE TOKEN AND INSERT RECORD
        token = secrets.token_urlsafe()
        try:
            new_user = User(id=uid, jaws_token=token, email=email, is_admin=admin)
            db.session.add(new_user)
            db.session.commit()
            logger.info(f"Added new user {uid} ({email})")
            print(f"User's access token:\n{token}")
        except Exception as e:
            logger.exception(f"Failed to add user: {e}")
            raise e


def jaws():
    """Entrypoint for jaws-server app."""
    cli()
