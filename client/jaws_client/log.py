import logging


def setup_logger(name: str, log_file=None, log_level="INFO") -> logging:
    """Configure logger.

    :param name: Name of the logger (e.g. __package__)
    :type name: str
    :param log_file: Path to write log file
    :type log_file: str
    :param log_level: Level at which to output logs
    :type log_level: str
    :return: logging object
    :rtype: obj
    """
    valid_log_levels = ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL")
    if log_level not in valid_log_levels:
        raise ValueError(
            f"Invalid log level: {log_level}; valid levels are {valid_log_levels}"
        )
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(levelname)s - %(module)s - %(message)s"
    )
    handler_stderr = logging.StreamHandler()
    handler_stderr.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(log_level)
    logger.addHandler(handler_stderr)

    if log_file:
        from logging.handlers import RotatingFileHandler
        handler_file = RotatingFileHandler(log_file, maxBytes=1024, backupCount=1, mode="a")
        handler_file.setFormatter(formatter)
        logger.addHandler(handler_file)
    return logger
