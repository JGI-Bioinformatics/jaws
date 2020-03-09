import logging
import os
import configparser


class Singleton(type):

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class JawsConfig(metaclass=Singleton):
    """Configuration singleton class"""
    config = None
    db = None
    session = None

    def __init__(self, config_file=None):
        """Constructor

        :param config_file: Path to configuration file in INI format
        :type config_file: str
        """
        logger = logging.getLogger(__package__)
        logger.debug('loading configuration...')
        if not config_file:
            raise FileNotFoundError("config file not specified")
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"{config_file} does not exist")
        self.config = configparser.ConfigParser()
        self.config.read(config_file)

    def get(self, section, key, default=None):
        """Get a configuration value.

    :param section: top-level section of the config
    :type section: str
    :param key: second-level parameter key
    :type key: str
    :return: the value is always a string; typecast as necessary
    :rtype: str
    """
        return self.config.get(section, key, default)


conf = None
