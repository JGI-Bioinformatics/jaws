#!/usr/bin/env python3

"""
JAWS Central REST Server
"""

import os
import click
import connexion
from jaws_central import config, log, rpc_manager, database


@click.group()
def central() -> None:
    pass


@central.command()
@click.option("--config", "config_file", default="jaws_central.yml", help="Central configuration YAML file")
@click.option("--log", "log_file", default="jaws-central.log", help="Log file")
def serve(config_file: str, log_file: str) -> None:
    """
    Start JAWS-Central REST server.
    """
    logger = log.setup_logger(__package__, log_file)
    logger.debug("Starting jaws-central server")

    config.conf = config.JawsConfig(config_file)
    database.db = database.JawsDb(config.conf)
    rpc_manager.rpc = rpc_manager.JawsRpc(config.conf)

    logger.debug("Initializing Connexion app")
    basedir = os.path.abspath(os.path.dirname(__file__))
    connex = connexion.FlaskApp(__name__, specification_dir=basedir)
    try:
        connex.add_api("swagger.yml")
    except Exception as e:
        raise e("Unable to initialize Connexion app")
    connex.run(host="0.0.0.0", port=5000, debug=False)
