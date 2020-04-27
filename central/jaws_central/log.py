import logging


def setup_logger(name: str, log_file="jaws_central.log") -> logging:
    """Configure the logging singleton object for the package.

    :param log_file: Path to the log file
    :type log_file: str
    :return: logging object which is a singleton
    :rtype: obj
    """
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(levelname)s - %(module)s - %(message)s"
    )

    handler_stderr = logging.StreamHandler()
    handler_stderr.setFormatter(formatter)

    handler_file = logging.FileHandler(log_file)
    handler_file.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler_stderr)
    logger.addHandler(handler_file)
    return logger
