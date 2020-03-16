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

    def __init__(self, config_file=None):
        """Constructor

        :param config_file: Path to configuration file in INI format
        :type config_file: str
        """
        self.logger = logging.getLogger(__package__)
        self.logger.debug('loading configuration...')
        if not config_file:
            raise FileNotFoundError("config file not specified")
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"{config_file} does not exist")
        self.config = configparser.ConfigParser()
        self.config.read(config_file)

    def get(self, section, key, default=None):
        """Get a configuration value.

    :param section: name of config section
    :type section: str
    :param key: parameter key
    :type key: str
    :return: the value is always a string; typecast as necessary
    :rtype: str
    """
        if section not in self.config:
            self.logger.warn(f"Config file missing section {section}")
            return default
        if key not in self.config[section]:
            self.logger.warn(f"Config file missing param {section}/{key}")
            return default
        return self.config.get(section, key)


conf = None
