#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Copyright 2018-2019 (C) DOE JGI, LBL
# Seung-Jin Sul (ssul@lbl.gov)

from __future__ import print_function
import sys
import os
import argparse
import getpass
import shutil
import fileinput

import uuid
import shortuuid
if sys.version_info[0] < 3:
    import cPickle
else:  # py3
    import _pickle as cPickle
import zlib

import multiprocessing
import threading
import functools
import time
import signal

import datetime
import socket
import subprocess
#import atexit

import json
import pprint

# import sqlite3 as sqlite
# import psutil
import logging

try:
    # Supports v0.12.0 and higher. The latest is v1.1.0 @ Aug 2019
    import pika
except ImportError:
    sys.stderr.write("Exception: Error importing pika. Please check your installation.")
    sys.exit(1)

from Utils.MsgCompress import zdumps, zloads
from Utils.Run import makeDir, run_sh_command, make_dir_p
from Config import *

logger = logging.getLogger(__name__)

# coloredlogs.install(level='DEBUG')
# coloredlogs.install(level='DEBUG', logger=logger)

## To hide "No handlers could be found for logger "pika.adapters.base_connection"" warning
logging.getLogger('pika').setLevel(logging.INFO)
# logging.getLogger('pika').setLevel(logging.DEBUG)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


#-------------------------------------------------------------------------------
def setup_custom_logger(level, logDestDir, bStreamLogging=True, bFileLogging=False, workerId=None):
#-------------------------------------------------------------------------------
    """
    Setting up logging

    @param level: logger level
    """
    ## Set custom loglevel name for printing resource usage
    ## CRITICAL = 50, ERROR = 40, WARNING = 30, INFO = 20, DEBUG = 10, NOTSET = 0.
    DEBUG_LEVELV_NUM = 60
    logging.addLevelName(DEBUG_LEVELV_NUM, "RESOURCE")

    def resource(self, message, *args, **kws):
        self._log(DEBUG_LEVELV_NUM, message, args, **kws)
    logging.Logger.resource = resource

    numericLevel = getattr(logging, level.upper(), None)
    if not isinstance(numericLevel, int):
        raise ValueError('Invalid log level: %s' % level)

    formatter = logging.Formatter('%(asctime)s | %(module)s | %(funcName)s | %(levelname)s : %(message)s')
    logger.setLevel(numericLevel)

    ## StreamLogger
    if bStreamLogging:
        streamLogger = logging.StreamHandler()
        streamLogger.setLevel(numericLevel)
        streamLogger.setFormatter(formatter)
        logger.addHandler(streamLogger)

    ## FileLogger
    if bFileLogging:
        dateStr = datetime.datetime.now().strftime("%Y-%m-%d")
        if workerId:
            dateStr = datetime.datetime.now().strftime("%Y-%m-%d")
            dateStr = workerId + '_' + dateStr

        if logDestDir:
            logDir = logDestDir
        else:
            logDir = '%s/logs' % (os.getcwd())

        makeDir(logDir)

        logFileName = '%s/jtm_%s.log' % (logDir, dateStr)
        if workerId:
            logFileName = '%s/jtm_worker_%s.log' % (logDir, dateStr)

        fileLogger = logging.FileHandler(logFileName)
        fileLogger.setFormatter(formatter)
        fileLogger.setLevel(numericLevel)
        logger.addHandler(fileLogger)
        logger.info("Log file name: %s" % (logFileName))


# class ConnWrapper(object):
#     _conn = None
#     def __init__(self, conn):
#         self._conn = conn