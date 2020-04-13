#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)
#
import datetime
import logging
import os

from jaws_jtm.lib.run import make_dir

logger = logging.getLogger(__package__)

# To hide "No handlers could be found for logger "pika.adapters.base_connection"" warning
logging.getLogger('pika').setLevel(logging.INFO)


# -------------------------------------------------------------------------------
def setup_custom_logger(level, log_dest_dir, b_stream_loggin=True, b_file_loggin=False, worker_id=None):
    """
    Setting up logging

    @param level: logger level
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
        raise ValueError('Invalid log level: %s' % level)

    formatter = logging.Formatter('%(asctime)s | %(module)s | %(funcName)s | %(levelname)s : %(message)s')
    logger.setLevel(numeric_level)

    # StreamLogger
    if b_stream_loggin:
        streamLogger = logging.StreamHandler()
        streamLogger.setLevel(numeric_level)
        streamLogger.setFormatter(formatter)
        logger.addHandler(streamLogger)

    # FileLogger
    if b_file_loggin:
        datetime_str = datetime.datetime.now().strftime("%Y-%m-%d")
        if worker_id:
            datetime_str = datetime.datetime.now().strftime("%Y-%m-%d")
            datetime_str = worker_id + '_' + datetime_str

        assert log_dest_dir
        worker_log_dir_name = os.path.join(log_dest_dir, "worker")

        make_dir(worker_log_dir_name)
        try:
            os.chmod(worker_log_dir_name, 0o775)
        except OSError:
            logger.warning("Cannot change the permission of {} to 0775".format(worker_log_dir_name))

        log_file_name = '%s/jtm_%s.log' % (worker_log_dir_name, datetime_str)
        if worker_id:
            log_file_name = '%s/jtm_worker_%s.log' % (worker_log_dir_name, datetime_str)

        file_logger = logging.FileHandler(log_file_name)
        file_logger.setFormatter(formatter)
        file_logger.setLevel(numeric_level)
        logger.addHandler(file_logger)
        logger.info("Log file name: %s" % (log_file_name))
        try:
            os.chmod(log_file_name, 0o775)
        except OSError:
            logger.warning("Cannot change the permission of {} to 0775".format(log_file_name))
            raise
