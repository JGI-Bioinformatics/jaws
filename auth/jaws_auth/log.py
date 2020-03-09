import logging


def setup_logger(name, log_file="jaws_central.log"):
    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')

    handler_stderr = logging.StreamHandler()
    handler_stderr.setFormatter(formatter)

    handler_file = logging.FileHandler(log_file)
    handler_file.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler_stderr)
    logger.addHandler(handler_file)
    return logger
