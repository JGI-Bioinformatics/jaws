#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)
#
import logging
import os
from logging import handlers
from jaws_jtm.lib.run import make_dir


logger = logging.getLogger(__package__)

# To hide "No handlers could be found for logger "pika.adapters.base_connection"" warning
logging.getLogger("pika").setLevel(logging.INFO)


# -------------------------------------------------------------------------------
def setup_custom_logger(
    level,
    log_dest_dir,
    log_file_name,
    b_stream_logging=True,
    b_file_logging=False,
    worker_id=None,
):
    """
    Setting up logging

    @param level: logger level
    @param log_dest_dir: log dir path
    @param log_file_name: log file name (full path)
    @param b_stream_logging: enable stream logging
    @param b_file_logging: enable file logging
    @param worker_id: if specidied, use worker_id to create a log file name
    """
    # Set custom loglevel name for printing resource usage
    # CRITICAL = 50, ERROR = 40, WARNING = 30, INFO = 20, DEBUG = 10, NOTSET = 0.
    DEBUG_LEVELV_NUM = 60
    logging.addLevelName(DEBUG_LEVELV_NUM, "RESOURCE")

    def resource(self, message, *args, **kws):
        self._log(DEBUG_LEVELV_NUM, message, args, **kws)

    logging.Logger.resource = resource

    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % level)

    formatter = logging.Formatter(
        "%(asctime)s | %(module)s | %(lineno)d | %(funcName)s | %(levelname)s : %(message)s",
        "%Y-%m-%d %H:%M:%S",
    )
    logger.setLevel(numeric_level)

    # StreamLogger
    if b_stream_logging:
        streamLogger = logging.StreamHandler()
        streamLogger.setLevel(numeric_level)
        streamLogger.setFormatter(formatter)
        logger.addHandler(streamLogger)

    # FileLogger
    if b_file_logging:
        assert log_dest_dir
        make_dir(log_dest_dir)
        try:
            os.chmod(log_dest_dir, 0o775)
        except OSError:
            logger.warning(
                "Cannot change the permission of {} to 0775".format(log_dest_dir)
            )

        # Rotational log: 100MB each, total 4 log files
        file_logger = handlers.RotatingFileHandler(
            log_file_name, maxBytes=100000000, backupCount=3
        )
        file_logger.setFormatter(formatter)
        file_logger.setLevel(numeric_level)
        logger.addHandler(file_logger)
        logger.info("Log file name: %s" % (log_file_name))
        try:
            os.chmod(log_file_name, 0o775)
        except OSError:
            logger.warning(
                "Cannot change the permission of {} to 0775".format(log_file_name)
            )
            raise
