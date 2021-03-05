#!/usr/bin/env python
"""
FLASK application containing endpoints to monitor site specific entities. The following endpoints are
available:

- Reports the percent of free disk space for the given <file_path>
  Ex: localhost:<port>/disk_free_pct/<file_path>
  Returns dictionary:
  {
      "disk_free_pct": float value
  }

- Reports the pid of a supervisord
  Ex: localhost:<port>/supervisor_pod?config=name of supervisor config&cmd=name and path to supervisord
  Returns pid as int

- Reports the status of each process managed by a supervisord; status=1 or 0.
  localhost:<port>/supervisor_pod?config=name of supervisor config&cmd=name and path to supervisorctl
  Returns dictionary where key=name of supervisor process, value=status
"""

import os
import argparse
import shutil
import subprocess
import getpass
import re
from flask import Flask, request

app = Flask(__name__.split('.')[0])


def check_path(file_path):
    """Check if input file path is valid. Returns True or False.

    :param none
    :type none
    :return: argparser object
    :rtype: object
    """

    return True if os.path.isdir(file_path) else False


def run_command(cmd):
    """Execute input command, returning stdout, stderr, exitcode.

    :param cmd: command to execute
    :type cmd: str
    :return: stdout, stderr, exitcode
    :rtype: str, str, int
    """

    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    exitcode = process.returncode

    return stdout.decode('utf-8').strip(), stderr.decode('utf-8').strip(), exitcode


def rename_supervisor_process(name):
    """Given a process name, look for name that contains a colon. If found, remove the word
    before the colon along with the colon and return that name. If no colon is found
    in the name, return the original name.

    EX: name=jaws-dev:jaws-central-auth, return jaws-central-auth
        name=jaws-central-auth, return jaws-central-auth

    :param name: process name from supervisorctl
    :type name: string
    :return name: modified name
    :rtype name: string
    """

    rex = re.search('.+:(.+)', name)
    if rex:
        return rex.group(1)
    else:
        return name


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
    if not check_path(file_path):
        file_path = os.path.join('/', file_path)

    # if path exists, compute pct disk free
    if check_path(file_path):
        total, _, free = shutil.disk_usage(file_path)
        disk_free_pct = "%.1f" % (free/total*100)

    return {'disk_free_pct': float(disk_free_pct)}, 200


@app.route('/supervisor_pid')
def get_supervisor_pid():
    """Check if supervisord is running and return the pid if found.

    The flask endpoint for this url request is
    localhost:port/supervisor_pid?config=supervisor.conf&/path/supervisor/venv

    If config parameter contains a path, the basename is used.

    :param config: name of supervisor config file
    :type config: string
    :param cmd: path and name to supervisord cmd
    :type cmd: string
    :return: dictionary containing key='pid', value=pid
    :rtype: dictionary
    """

    pid = 0
    username = getpass.getuser()
    config = request.args.get('config', None)   # supervisord config file
    supervisor_cmd = request.args.get('cmd', None)  # supervisord cmd full path and name

    if not config:
        return {'pid': 0, 'error': 'Missing config parameter in url.'}, 400
    if not supervisor_cmd:
        return {'pid': 0, 'error': 'Missing cmd parameter in url.'}, 400
    if not os.path.isfile(supervisor_cmd):
        return {'pid': 0, 'error': f'Cannot find supervisor cmd: {supervisor_cmd}'}, 400

    config = os.path.basename(config)
    cmd = "ps aux | grep {0} | grep {1} | grep {2} | grep -v grep".format(
        username, supervisor_cmd, config)

    stdout, stderr, exitcode = run_command(cmd)

    if exitcode:
        return {'pid': 0, 'error': 'ps encountered an error'}, 400
    else:
        # parse pid from ps command. output looks like:
        # strong   48047  0.0  0.0  46548 22632
        values = stdout.split()
        if len(values) > 1:
            pid = int(values[1])

    return {'pid': pid}, 200


@app.route('/supervisor_process')
def get_supervisor_processes():
    """Lookup all processes managed by a supervisord and report the status of each process.
    Returns a dictionary where the key is the name of the supervisor process, the value is either 1 if
    up or 2 if down.

    :param config: name of supervisor config file
    :type config: string
    :param cmd: name and path to supervisorctl
    :type cmd: string
    :return: dictionary containing key=name of supervisor process, value=status as 0|1
    :rtype: dictionary
    """

    processes = {}
    config = request.args.get('config', None)  # supervisord config file
    supervisor_cmd = request.args.get('cmd', None)  # supervisord cmd full path and name

    if not config:
        return {'pid': 0, 'error': 'Missing config parameter in url.'}, 400
    if not supervisor_cmd:
        return {'pid': 0, 'error': 'Missing cmd parameter in url.'}, 400
    if not os.path.isfile(supervisor_cmd):
        return {'pid': 0, 'error': f'Cannot find supervisor cmd: {supervisor_cmd}'}, 400

    cmd = f"{supervisor_cmd} -c {config} status"

    stdout, _, exitcode = run_command(cmd)

    # Output of cmd should look something like:
    # jaws-dev:jaws-central-auth       RUNNING   pid 43761, uptime 0:09:17
    # jaws-dev:jaws-central-rest       RUNNING   pid 43759, uptime 0:09:17
    # jaws-dev:jaws-central-rpc        RUNNING   pid 43760, uptime 0:09:17
    # jaws-dev:jaws-site-central-rpc   RUNNING   pid 43758, uptime 0:09:17
    # jaws-dev:jaws-site-daemon        RUNNING   pid 43756, uptime 0:09:17
    # jaws-dev:jaws-site-jtm-rpc       RUNNING   pid 43757, uptime 0:09:17

    # supervisorctl exits with code 4 if error is 'http://localhost:nnnn refused connection'
    if exitcode == 4:
        return {'error': stdout}, 400

    # Parse output getting the name of the process minus any word before the colon and
    # the status. If the status is RUNNING, set status=1, else 0
    #
    if stdout:
        for line in stdout.split('\n'):
            line = line.strip()
            rows = line.split()
            if len(rows) > 1:
                name = rename_supervisor_process(rows[0])
                status = 1 if rows[1] == 'RUNNING' else 0
                processes[name] = status

    return processes, 200


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
