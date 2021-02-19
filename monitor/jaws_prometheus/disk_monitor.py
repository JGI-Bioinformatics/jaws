#!/usr/bin/env python

import os
import sys
import shutil
import connexion


def disk_free(filepath='/global/cscratch1'):
    total, use, free = shutil.disk_usage(filepath)
    return float("%.1f"%(free/total*100))


def main():
    basedir = os.path.abspath(os.path.dirname(__file__))
    port = 5000
    host = '0.0.0.0'
    app = connexion.App(__name__, specification_dir=basedir)
    app.add_api("disk_monitor_swagger.yml")
    app.run(host=host, port=port, debug=False)

if __name__ == '__main__':
    main()