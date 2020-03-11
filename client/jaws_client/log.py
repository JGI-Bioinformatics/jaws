import logging
from logging.handlers import RotatingFileHandler


def setup_logger(name: str, log_file=None) -> logging:
    """Configure logger.

    :param name: Name of the logger (e.g. __package__)
    :type name: str
    :param log_file: Path to write log file
    :type log_file: str
    :return: logging object
    :rtype: obj
    """
    if log_file is None:
        log_file = f"{name}.log"
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(levelname)s - %(module)s - %(message)s"
    )

    handler_stderr = logging.StreamHandler()
    handler_stderr.setFormatter(formatter)

    handler_file = RotatingFileHandler(log_file, maxBytes=512, backupCount=2, mode="a")
    handler_file.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler_stderr)
    logger.addHandler(handler_file)

    return logger
