#!/usr/bin/env python

"""
Miscellaneous utilities
"""

import sys
import os
import configparser
from flask import make_response, abort, Flask, request, redirect, url_for, current_app
import rpc_client


DEBUG = True if "JAWS_DEBUG" in os.environ else False

# JAWS-SITE CONFIG
if "JAWS_SITES_CONFIG" not in os.environ: sys.exit('Env var $JAWS_SITES_CONFIG not defined')
site_config = configparser.ConfigParser()
site_config.read_file(open(os.environ["JAWS_SITES_CONFIG"]))





def status(user):
    """
    Check system health
    """
    #config = current_app.config
    result = {
        "JAWS-Central": "UP",
        "JAWS-Catalog" : "UP",
        "JAWS-Auth" : "UP"
    }

    # INIT SITE RPC OBJECTS AND CHECK STATUS OF EACH
    for site_name in site_config.sections():
        if DEBUG: print("Initializing RPC client for %s" % (site_name,))
        rpc = rpc_client.RPC_Client(
            site_config[site_name]["amqp_host"],
            site_config[site_name]["amqp_user"],
            site_config[site_name]["amqp_password"],
            site_config[site_name]["amqp_queue"],
            vhost=site_config[site_name]["amqp_vhost"] )
        response = rpc.request("server_status")
        if "error" not in response:
            result[site_name+"-Site"] = "UP"
            result[site_name+"-Cromwell"] = "UP"
        elif response["error"]["code"] == 500:
            result[site_name+"-Site"] = "DOWN"
            result[site_name+"-Cromwell"] = "Unknown"
        else:
            result[site_name+"-Site"] = "UP"
            result[site_name+"-Cromwell"] = "DOWN"

    return result, 200

