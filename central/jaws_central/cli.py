#!/usr/bin/env python3

"""JAWS Central CLI"""

import os
import click
import connexion
from urllib.parse import quote_plus
import secrets
from jaws_central import config, log, rpc_manager
from jaws_central.models import db


@click.group()
def cli() -> None:
    pass


def find_config_file(env_var: str, config_file: str) -> str:
    """Find config file by checking env var, cwd, and home, in that order.

    :param config_file: filename of configuration file
    :type config_file: str
    :param env_var: name of config file environment variable
    :type env_var: str
    :return: absolute path if found, None otherwise.
    :rtype: str
    """
    if env_var in os.environ:
        cf = os.environ[env_var]
        if os.path.isfile(cf):
            return cf
        else:
            raise IOError(f"File not found: {cf}")

    cf = os.path.abspath(config_file)
    if os.path.isfile(cf):
        return cf

    cf = os.path.join(os.environ["HOME"], config_file)
    if os.path.isfile(cf):
        return cf
    raise IOError("Unable to find config file")


@cli.command()
@click.option(
    "--config", "config_file", default="jaws-central.ini", help="Config INI file"
)
@click.option("--log", "log_file", default="jaws-auth.log", help="Log file")
def auth(config_file, log_file):
    """Start JAWS OAuth server"""
    logger = log.setup_logger(__package__, log_file)
    logger.debug("Starting jaws-auth server")

    config_file = find_config_file("JAWS_AUTH_CONFIG", config_file)
    conf = config.JawsConfig(config_file)

    logger.debug("Initializing Connexion app")
    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp(__name__, specification_dir=basedir)
    connex.add_api("swagger.auth.yml")

    connex.app.config["SQLALCHEMY_DATABASE_URI"] = "%s://%s:%s@%s/%s" % (
        conf.get("DB", "dialect"),
        conf.get("DB", "user"),
        quote_plus(conf.get("DB", "password")),
        conf.get("DB", "host"),
        conf.get("DB", "db"),
    )
    connex.app.config["SQLALCHEMY_ECHO"] = False
    connex.app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False
    connex.app.config["SQLALCHEMY_ENGINE_OPTIONS"] = {
        "pool_pre_ping": True,
        "pool_recycle": 300,
    }
    db.init_app(connex.app)

    connex.run(host="0.0.0.0", port=3000, debug=False)


@cli.command()
@click.option(
    "--config",
    "config_file",
    default="jaws-central.ini",
    help="Central configuration INI file",
)
@click.option("--log", "log_file", default="jaws-rest.log", help="Log file")
def rest(config_file: str, log_file: str) -> None:
    """Start JAWS REST server."""
    logger = log.setup_logger(__package__, log_file)
    logger.debug("Starting jaws-central server")

    config_file = find_config_file("JAWS_REST_CONFIG", config_file)
    conf = config.JawsConfig(config_file)

    logger.debug("Initializing Connexion app")
    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp("JAWS_REST", specification_dir=basedir)
    connex.add_api("swagger.rest.yml")

    connex.app.config["SQLALCHEMY_DATABASE_URI"] = "%s://%s:%s@%s/%s" % (
        conf.get("DB", "dialect"),
        conf.get("DB", "user"),
        quote_plus(conf.get("DB", "password")),
        conf.get("DB", "host"),
        conf.get("DB", "db"),
    )
    connex.app.config["SQLALCHEMY_ECHO"] = False
    connex.app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False
    connex.app.config["SQLALCHEMY_ENGINE_OPTIONS"] = {
        "pool_pre_ping": True,
        "pool_recycle": 300,
    }
    db.init_app(connex.app)

    rpc_manager.rpc = rpc_manager.JawsRpc(conf)

    connex.run(host="0.0.0.0", port=5000, debug=False)


@cli.command()
@click.option(
    "--config", "config_file", default="jaws-central.ini", help="Config INI file"
)
@click.option("--log", "log_file", default="jaws-auth.log", help="Log file")
def create_tables(config_file: str, log_file: str) -> None:
    """Create all database tables"""
    logger = log.setup_logger(__package__, log_file)

    config_file = find_config_file("JAWS_REST_CONFIG", config_file)
    conf = config.JawsConfig(config_file)

    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp("JAWS_REST", specification_dir=basedir)
    connex.add_api("swagger.rest.yml")

    connex.app.config["SQLALCHEMY_DATABASE_URI"] = "%s://%s:%s@%s/%s" % (
        conf.get("DB", "dialect"),
        conf.get("DB", "user"),
        quote_plus(conf.get("DB", "password")),
        conf.get("DB", "host"),
        conf.get("DB", "db"),
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
@click.option(
    "--config", "config_file", default="jaws-central.ini", help="Config INI file"
)
@click.option("--log", "log_file", default="jaws-auth.log", help="Log file")
def add_user(
    uid: str, email: str, config_file: str, log_file: str, admin: bool = False
) -> None:
    """Add user and generate OAuth2 token."""
    logger = log.setup_logger(__package__, log_file)

    config_file = find_config_file("JAWS_REST_CONFIG", config_file)
    conf = config.JawsConfig(config_file)

    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp("JAWS_REST", specification_dir=basedir)
    connex.add_api("swagger.rest.yml")

    connex.app.config["SQLALCHEMY_DATABASE_URI"] = "%s://%s:%s@%s/%s" % (
        conf.get("DB", "dialect"),
        conf.get("DB", "user"),
        quote_plus(conf.get("DB", "password")),
        conf.get("DB", "host"),
        conf.get("DB", "db"),
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
