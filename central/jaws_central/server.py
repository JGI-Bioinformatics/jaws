#!/usr/bin/env python3

"""
JAWS Central REST Server
"""

import logging
import config

logging.basicConfig(level=logging.INFO)

connex = config.connex
connex.add_api('swagger.yml')

if __name__ == "__main__":
    connex.run(host='0.0.0.0', port=5000, debug=False)
