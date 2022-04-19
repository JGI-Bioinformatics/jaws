#!/usr/bin/env python
"""
Condor pool manager

Author: Seung-Jin Sul (ssul@lbl.gov)

Steps
1. Collect the IDLE condor job IDs
2. Get the number of nodes required per memory requirement (normal, jgi_shared, jgi_exvivo)
   normal: <= 118G
   jgi_shared: 118 ~ 740G
   jgi_exvivo: 740 ~ 1450G
3. Execute sbatch if the current SLURM jobs + number of sbatches needed < MAX pool size
4. Check if MIN pool size is being maintained, if not add more nodes up to MIN

"""
import click
import logging
import os
import time
import shlex
from datetime import datetime
from jaws_condor import config
from jaws_condor.utils import run_sh_command, run_slurm_cmd


logger = None


def start_file_logger(filename, name='jaws_condor_pool_add.log', level=logging.DEBUG, format_string=None):
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


def add():
    pass
