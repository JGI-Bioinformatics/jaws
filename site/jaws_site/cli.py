#!/usr/bin/env python

"""
JAWS Site server runs at each computing site and is comprised of:
(a) RPC server for handling user (sync) requests and
(b) daemon for performing periodic (async) maintenance tasks.
Each computing site also has a Cromwell server instance, typically installed on the same server.
"""

import os
import sys
import click
from jaws_site import config, database, rpc_server, jawsd, log, dispatch


log_file = os.environ["JAWS_SITE_LOG"] if "JAWS_SITE_LOG" in os.environ else "./jaws-site.log"
logger = log.setup_logger(__package__, log_file)


def find_config_file(env_var, config_file):
    """Find config file by checking env var, cwd, and home, in that order.

    :param config_file: filename of configuration file
    :type config_file: str
    :param env_var: name of config file environment variable
    :type env_var: str
    :return: absolute path if found, None otherwise.
    :rtype: str
    """
    if env_var is os.environ:
        cf = os.environ[env_var]
        if os.path.isfile(cf):
            return cf

    cf = os.path.abspath(config_file)
    if os.path.isfile(cf):
        return cf

    cf = os.path.join(os.environ["HOME"], config_file)
    if os.path.isfile(cf):
        return cf
    return None


config_file = find_config_file("JAWS_SITE_CONFIG", "jaws-site.ini")
if config_file is None:
    sys.exit("Configuration required")
conf = config.JawsConfig(config_file)


db = database.JawsDb(conf)


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def serve():
    """JAWS-Site"""
    pass


@serve.command()
def server():
    """Start JAWS-Site RPC-server."""
    app = rpc_server.RpcServer(conf)
    app.start_server()


@serve.command()
def daemon():
    """Start JAWS-Site daemon."""
    daemon = jawsd.JAWSd(conf, db)
    daemon.start_daemon()


cli = click.CommandCollection(sources=[serve, dispatch.cli, jawsd.cli])

if __name__ == "__main__":
    cli()
