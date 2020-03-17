import logging
import os
import configparser


conf = None


class Singleton(type):

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class ConfigurationError(Exception):
    pass


class JawsConfig(metaclass=Singleton):
    """Configuration singleton class"""
    config: configparser

    def __init__(self, config_file: str = None):
        """Constructor sets global singleton.

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
        for section in ['AMQP', 'RPC', 'GLOBUS', 'DB', 'CROMWELL']:
            if section not in self.config:
                raise ConfigurationError(f"Config missing required section: {section}")
        global conf
        conf = self

    def get(self, section: str, key: str) -> str:
        """Get a configuration value.

        :param section: name of config section
        :type section: str
        :param key: parameter key
        :type key: str
        :return: the value is always a string; typecast as necessary
        :rtype: str
        """
        if section in self.config:
            if key in self.config[section]:
                return self.config[section][key]
            else:
                raise ConfigurationError(f"Section {section} does not have {key} parameter")
        else:
            raise ConfigurationError(f"Config file does not have {section} section")
