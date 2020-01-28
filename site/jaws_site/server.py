#!/usr/bin/env python

from jaws_site import rpc_server, config

config = config.Config(env="JAWS_SITE_CONFIG")

if __name__ == '__main__':
    print("Starting RPC server")
    app = rpc_server.RpcServer(config)
    app.start_server()
