#!/usr/bin/env python

"""
JAWS Site server runs at each computing site and is comprised of:
(a) RPC server for handling user (sync) requests and
(b) daemon for performing periodic (async) maintenance tasks.
Each computing site also has a Cromwell server instance, typically installed on the same server.
"""

import os
import click
from jaws_site import database, rpc_server, jawsd, log
from jaws_site.config import jaws_config


@click.group()
def site():
    pass


@site.command()
@click.option("--config", "config_file", default="jaws-site.ini", help="Config INI file")
@click.option("--log", "log_file", default="jaws-site.log", help="Log file")
def serve(log_file):
    """
    Start RPC-server and daemon.
    """
    logger = log.setup_logger(__package__, log_file)
    logger.debug("Starting jaws-site server")
    db = database.JawsDb(jaws_config)
    if os.fork():
        app = rpc_server.RpcServer()
        app.start_server()
    else:
        jd = jawsd.JAWSd(db)
        jd.start_daemon()
