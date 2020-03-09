#!/usr/bin/env python3

"""
JAWS Auth Server
"""

import os
import click
import connexion
from jaws_auth import config, log


@click.group()
def auth():
    pass


@auth.command()
@click.option("--config", "config_file", default="jaws-auth.ini", help="Config INI file")
@click.option("--log", "log_file", default="jaws-auth.log", help="Log file")
def serve(config_file, log_file):
    """
    Start OAuth2 server
    """
    logger = log.setup_logger(__package__, log_file)
    logger.debug("Starting jaws-auth server")
    conf = config.JawsConfig(config_file)
    conf.init_db()

    logger.debug("Initializing Connexion app")
    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp(__name__, specification_dir=basedir)
    try:
        connex.add_api("swagger.yml")
    except Exception as e:
        raise e("Unable to initialize Connexion app")
    connex.run(host="0.0.0.0", port=3000, debug=False)
