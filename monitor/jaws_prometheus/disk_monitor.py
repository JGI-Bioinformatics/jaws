#!/usr/bin/env python
"""
FLASK application to report the percent of free disk space for a given
file directory. The flask endpoint for checking is
localhost:<port>/disk_free_pct/<file_path>

Example:
localhost:51000/disk_free_pct/global/cscratch

Returns a json data in the format:
{
    "disk_free_pct": float value
}
"""

import os
import argparse
import shutil

from flask import Flask
app = Flask(__name__)


def _check_path(file_path):
    """Check if input file path is valid. Returns True or False.

    :param none
    :type none
    :return: argparser object
    :rtype: object
    """
    return True if os.path.isdir(file_path) else False


@app.route('/disk_free_pct/<path:file_path>')
def disk_free_pct(file_path):
    """Check the percent of free disk for the given path and return value as a float.
    If input file path is invalid, returns -1.

    The flask endpoint for this url request is
    localhost:port/disk_free_pct/<file_path>

    EX: localhost:51000/disk_free_pct/global/cscratch

    :param file_path: file path to check for free disk percent.
    :type file_path: string
    :return: dictionary containing key='disk_free_pct', value=float
    :rtype: dictionary
    """
    disk_free_pct = -1

    # if path doesn't exists, try prepend path with '/'
    if not _check_path(file_path):
        file_path = os.path.join('/', file_path)

    # if path exists, compute pct disk free
    if _check_path(file_path):
        total, use, free = shutil.disk_usage(file_path)
        disk_free_pct = "%.1f" % (free/total*100)

    return {'disk_free_pct': float(disk_free_pct)}


def get_args():
    """Parse command line arguments.

    :param none
    :type none
    :return: argparser object
    :rtype: object
    """
    prog_desc = '''
    '''
    parser = argparse.ArgumentParser(description=prog_desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--host', help='host name', default='0.0.0.0')
    parser.add_argument('--port', help='port', type=int, default=51000)
    parser.add_argument('--debug', action='store_true', help='debug')
    return parser.parse_args()


def main():
    """Main routine.
    """
    args = get_args()
    app.run(args.host, args.port, args.debug)


if __name__ == '__main__':
    main()
