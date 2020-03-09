"""
Configuration singleton, loads values from provided INI infile.
"""

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
    config = None
    db = None
    session = None

    def __init__(self, config_file=None):
        """Constructor

        :param config_file: path to configuration file in INI format
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

    def get(self, section, key):
        """Get configuration value

        :param section: top-level section of config file
        :type section: str
        :param key: second-level section of config file
        :type key: str
        :return: The return value is always a string; typecast as necessary
        :rtype: str
        """
        return self.config.get(section, key)


conf = None
