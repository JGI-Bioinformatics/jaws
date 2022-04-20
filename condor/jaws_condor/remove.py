#!/usr/bin/env python
"""
Condor pool manager

Author: Seung-Jin Sul (ssul@lbl.gov)

Steps
1. Collect all RUNNING and PENDING SLURM job ids for Condor pool
2. Get each number of RUNNING and PENDING jobs per the memory requirement (normal, jgi_shared, jgi_exvivo)
3. Check if any of running SLURM jobs are using the node(s) and `scancel` them if not and the number of R and PD SLURM
   jobs is bigger than MIN pool size

"""
import click
import logging
import os
import time
import shlex
from datetime import datetime
import configparser
import jaws_condor.config
from jaws_condor import config
form jaws_condor.utils import run_sh_command


logger = None


def start_file_logger(filename, name='jaws_condor_pool_remove.log', level=logging.DEBUG, format_string=None):
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


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default=None,
                        help="Config file")
    parser.add_argument("-l", "--logfile", default=None,
                        help="Path to logfile")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")

    args = parser.parse_args()
    jaws_condor.config.Configuration(args.config)
    site_id = jaws_condor.config.conf.get_site_id()
    site_config = configparser.ConfigParser()
    if site_id == "CORI":
        site_config.read("condor_config/cori_config.ini")
    elif site_id == "JGI":
        site_config.read("condor_config/jgi_config.ini")
    elif site_id == "TAHOMA":
        site_config.read("condor_config/tahoma_config.ini")
    else:
        raise ValueError("Unknown side_id specified in condor backend config")
