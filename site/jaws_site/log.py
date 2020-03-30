import logging


def setup_logger(name, log_file=None):
    """Configure logging singleton for the package.

    :param log_file: Path to log output file
    :type log_file: str
    :return: logging object, a singleton
    :rtype: obj
    """
    if log_file is None:
        log_file = f'{name}.log'
    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')

    handler_stderr = logging.StreamHandler()
    handler_stderr.setFormatter(formatter)

    handler_file = logging.FileHandler(log_file)
    handler_file.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.addHandler(handler_stderr)
    logger.addHandler(handler_file)

    return logger
