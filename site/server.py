#!/usr/bin/env python

"""
JAWS Site server runs at each computing site and is comprised of (a) RPC server for handling user requests and (b) daemon for performing periodic maintenance tasks.
Each computing site also has a Cromwell server instance, typically installed on the same server as JAWS-Site; communication between these services is via REST.

"""

import sys
import os
import schedule
import time
import configparser
from jaws_site import rpc_server, jawsd

# CONFIG FILE
if "JAWS_SITE_CONFIG" not in os.environ: sys.exit("Env var \$JAWS_SITE_CONFIG not defined")
config_file = os.environ["JAWS_SITE_CONFIG"]
config = configparser.ConfigParser()
config.read(config_file)

if __name__ == '__main__':
    """
    Fork and start RPC server (for responding to user requests) and daemon (for maintenance tasks).
    """
    if os.fork():
        print("Starting RPC server")
        app = rpc_server.RpcServer(config)
        app.start_server()
    else:
        print("Starting JAWS daemon")
        schedule.every(1).minutes.do(jawsd.check_uploads)
        schedule.every(1).minutes.do(jawsd.check_cromwell)
        while True:
            schedule.run_pending()
            time.sleep(10)
