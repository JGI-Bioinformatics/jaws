"""
Configuration singleton, loads values from provided INI infile.
"""

import logging
import os
import configparser


conf = None


class Singleton(type):

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super().__call__(*args, **kwargs)
        return cls._instances[cls]

    def _destructor(cls):
        if cls in cls._instances:
            del cls._instances[cls]


class ConfigurationError(Exception):
    def __init__(self, message):
        super().__init__(message)


class Configuration(metaclass=Singleton):
    """Configuration singleton."""

    defaults = {
        "USER": {
            "token": "",
            "staging_dir": ""
        },
        "JAWS": {
            "name": "JAWS",
            "site_id": "",
            "url": "http://localhost:5000",
            "womtool": ""
        },
        "GLOBUS": {
            "client_id": "",
            "endpoint_id": "",
            "basedir": "/"
        }
    }

    config = None

    def __init__(self, config_file) -> None:
        """Initialize the configuration object singleton

        :param config_file: Path to configuration file
        :type config_file: str
        """
        logger = logging.getLogger(__package__)
        logger.debug(f"Loading config from {config_file}")
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"{config_file} does not exist")
        self.config = configparser.ConfigParser()
        self.config.read_dict(self.defaults)
        try:
            self.config.read(config_file)
        except Exception as error:
            logger.exception(f"Unable to load config from {config_file}: {error}")
            raise
        global conf
        conf = self

    def get(self, section: str, key: str) -> str:
        """Get a configuration value (which is always a string).

        :param section: The section of the INI file
        :type section: str
        :param key: The parameter name
        :type key: str
        :return: Returns value, if exists; raises ConfigurationError otherwise.
        :rtype: str
        """
        if section not in self.config:
            raise ConfigurationError(
                f"Config file doesn't have section {section} defined"
            )
        elif key not in self.config[section]:
            raise ConfigurationError(
                f"Config file doesn't have parameter {section}/{key} defined"
            )
        return self.config[section][key]
