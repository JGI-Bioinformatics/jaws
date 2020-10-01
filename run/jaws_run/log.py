import logging


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
        log_file = f'{name}.log'
    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')

    handler_stderr = logging.StreamHandler()
    handler_stderr.setFormatter(formatter)

    handler_file = logging.FileHandler(log_file)
    handler_file.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(log_level)
    logger.addHandler(handler_stderr)
    logger.addHandler(handler_file)

    return logger
