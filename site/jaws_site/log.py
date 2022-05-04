import logging
from logging.handlers import RotatingFileHandler


def setup_logger(name, log_file=None, log_level="INFO"):
    """Configure logging singleton for the package.

    :param log_file: Path to log output file
    :type log_file: str
    :param log_level: Level at which to output logs
    :type log_level: str
    :return: logging object, a singleton
    :rtype: obj
    """
    valid_log_levels = ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL")
    if log_level not in valid_log_levels:
        raise ValueError(
            f"Invalid log level: {log_level}; valid levels are {valid_log_levels}"
        )

    if log_file is None:
        log_file = f"{name}.log"
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(levelname)s - %(module)s - %(message)s"
    )

    handler_stderr = logging.StreamHandler()
    handler_stderr.setFormatter(formatter)

    # Rotational log: 100MB each, total 4 log files
    handler_file = RotatingFileHandler(log_file, maxBytes=100000000, backupCount=3)
    handler_file.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(log_level)
    logger.addHandler(handler_stderr)
    logger.addHandler(handler_file)

    return logger
