#!/usr/bin/env python

"""jaws-auth : OAUTH service"""

import os
import click
import logging
import connexion
from urllib.parse import quote_plus
import secrets
from sqlalchemy.pool import QueuePool
from jaws_auth import config, log
from jaws_auth.models import db


JAWS_LOG_ENV = "JAWS_AUTH_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), f"{__package__}.log")
JAWS_CONFIG_ENV = "JAWS_AUTH_CONFIG"
JAWS_CWD_CONFIG = os.path.join(os.getcwd(), f"{__package__}.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(config_file: str, log_file: str, log_level: str) -> None:
    """JAWS-Auth"""
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
    """Start JAWS OAuth server"""
    logger = logging.getLogger(__package__)
    logger.debug("Initializing OAuth server")
    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp(__name__, specification_dir=basedir)
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
        "poolclass": QueuePool,
        "pool_pre_ping": True,
        "pool_recycle": 3600,
        "pool_size": 5,
        "max_overflow": 10,
        "pool_timeout": 30
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
    port = int(config.conf.get("AUTH", "port"))  # defaults to 3000

    # start OAuth server
    connex.run(host="0.0.0.0", port=port, debug=False)


@cli.command()
@click.argument("uid")
@click.option("--admin", is_flag=True, default=False, help="Grant admin privileges")
def add_user(
    uid: str, admin: bool = False
) -> None:
    """Add user and generate OAuth2 token."""
    logger = logging.getLogger(__package__)
    logger.debug(f"Adding new user, {uid}")

    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp(__name__, specification_dir=basedir)
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

    with connex.app.app_context():
        # CHECK IF UID EXISTS
        from jaws_auth.models import User
        user = db.session.query(User).get(uid)
        if user is not None:
            msg = f"Cannot add user {uid}; user.id already taken."
            logger.debug(msg)
            raise ValueError(msg)

        # GENERATE TOKEN AND INSERT RECORD
        token = secrets.token_urlsafe()
        try:
            new_user = User(id=uid, access_token=token, is_admin=admin)
            db.session.add(new_user)
            db.session.commit()
            logger.info(f"Added new user {uid} (is_admin={admin})")
            print(f"User's access token:\n{token}")
        except Exception as e:
            logger.exception(f"Failed to add user: {e}")
            raise e


def jaws():
    """Entrypoint for jaws-auth app."""
    cli()
