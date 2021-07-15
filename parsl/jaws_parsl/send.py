# run like:
# python send.py -c <# of CPUs> -m <memory> -s <site> -cmd $(pwd)/hello.sh

import click
import logging
import os
from jaws_parsl import config
from multiprocessing.connection import Client

logger = None


def start_file_logger(filename, name='jaws-parsl-send.log', level=logging.DEBUG, format_string=None):
    """Add a stream log handler.

    Args:
        - filename (string): Name of the file to write logs to
        - name (string): Logger name
        - level (logging.LEVEL): Set the logging level.
        - format_string (string): Set the format string

    Returns:
       -  None
    """
    if format_string is None:
        format_string = "%(asctime)s.%(msecs)03d %(name)s:%(lineno)d [%(levelname)s]  %(message)s"

    global logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename)
    handler.setLevel(level)
    formatter = logging.Formatter(format_string, datefmt='%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)


@click.command()
@click.option("-c", "--cpus", help="Number of CPUs needed")
@click.option("-m", "--memory", help="Amount of memory needed")
@click.option("-s", "--site", help="Compute site")
@click.option("-cmd", "--command", help="Command script")
@click.option("-l", "--loglevel", default="debug", help="Logging level ('debug' or 'info')")
@click.option("-f", "--logfile", default=None, help="Path to logfile")
def send(cpus, memory, site, command, loglevel, logfile):
    job_info = {
        "cpus": cpus,
        "memory": memory,
        "command": command,
        "site": site
    }

    if logfile:
        logfile_path = logfile
    else:
        logfile_path = './jaws-parsl-send.log'
    os.makedirs(os.path.dirname(logfile_path), exist_ok=True)
    start_file_logger(logfile_path, level=logging.DEBUG if loglevel == "debug" else logging.INFO)
    logger.info("Starting")

    mp_params = config.conf.get_mp_params()
    mp_password = mp_params["password"]
    mp_host = mp_params["host"]
    mp_port = mp_params["port"]

    conn = Client((mp_host, mp_port), authkey=bytes(mp_password, encoding='utf-8'))
    conn.send(job_info)
    conn.close()

    logger.debug(" [x] Sent job info.")


if __name__ == "__main__":
    send()
