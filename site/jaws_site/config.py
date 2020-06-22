import os
import logging
import configparser
import jaws_site.utils


conf = None


class ConfigurationError(Exception):
    def __init__(self, message):
        super().__init__(message)


class Configuration(metaclass=jaws_site.utils.Singleton):

    """Configuration singleton class"""

    defaults = {
        "SITE_RPC_SERVER": {
            "host": "localhost",
            "port": "5672",
            "vhost": "site",
            "user": "guest",  # default from docker container
            "password": "guest",  # default from docker container
            "num_threads": 5,
            "max_retries": 5},
        "CENTRAL_RPC_CLIENT": {
            "host": "localhost",
            "port": "5672",
            "vhost": "jaws_central",
            "user": "",
            "password": "",
        },
        "GLOBUS": {"client_id": "", "endpoint_id": "jaws-testing", "root_dir": "/"},
        "DB": {
            "host": "localhost",
            "port": "3306",
            "user": "jaws",
            "password": "jawstest",
            "db": "jaws-workflows",
            "dialect": "mysql+mysqlconnector",
        },
        "CROMWELL": {
            "url": "localhost:8000",
        },
        "SITE": {
            "id": "local",
            "staging_subdirectory": "staging",
            "results_subdirectory": "results",
        },
    }

    config = None

    def __init__(self, config_file: str) -> None:
        """Constructor sets global singleton.

        :param config_file: Path to configuration file in INI format
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
            logger.exception(f"Unable to load config file {config_file}: {error}")
            raise
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
        return self.config.get(section, key)

    def get_section(self, section: str) -> dict:
        """Get a configuration section.

        :param section: name of config section
        :type section: str
        :return: A copy of the requested section
        :rtype: dict
        """
        result = {}
        for key, value in self.config.items(section):
            result[key] = value
        return result
