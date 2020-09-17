#!/usr/bin/env python3

"""JAWS-User CLI"""

import os
import click
import logging
import connexion
from urllib.parse import quote_plus
import secrets
from jaws_user import config, log
from jaws_user.database import Session
from jaws_user.models import User
from jaws_rpc import rpc_index, rpc_server


JAWS_LOG_ENV = "JAWS_USER_LOG"
JAWS_CWD_LOG = os.path.join(os.getcwd(), f"{__package__}.log")
JAWS_CONFIG_ENV = "JAWS_USER_CONFIG"
JAWS_CWD_CONFIG = os.path.join(os.getcwd(), f"{__package__}.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(config_file: str, log_file: str, log_level: str) -> None:
    """JAWS-User"""
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
def server() -> None:
    """Start JAWS-User RPC server."""
    from jaws_user import rpc_operations

    rpc_params = config.conf.get_section("RPC_SERVER")
    app = rpc_server.RpcServer(rpc_params, rpc_operations.operations)
    app.start_server()


# @cli.command()
# @click.argument("uid")
# @click.argument("email")
# @click.option("--admin", is_flag=True, default=False, help="Grant admin privileges")
# def add_user(
#    uid: str, email: str, admin: bool = False
# ) -> None:
#    """Add user and generate OAuth2 token."""
#    logger = logging.getLogger(__package__)
#    logger.debug(f"Adding new user, {uid}")
#
#    # CHECK IF UID EXISTS
#    session = Session()
#    user = session.query(User).get(uid)
#    if user is not None:
#        msg = f"Cannot add user {uid}; user.id already taken."
#        logger.debug(msg)
#        raise ValueError(msg)
#
#    # GENERATE TOKEN AND INSERT RECORD
#    token = secrets.token_urlsafe()
#    try:
#        new_user = User(id=uid, jaws_token=token, email=email, is_admin=admin)
#        db.session.add(new_user)
#        db.session.commit()
#        logger.info(f"Added new user {uid} ({email})")
#        print(f"User's access token:\n{token}")
#    except Exception as e:
#        logger.exception(f"Failed to add user: {e}")
#        raise e


def jaws():
    """Entrypoint for jaws-user app."""
    cli()
