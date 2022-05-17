import os
import logging
import configparser
from jaws_site.env_interpolation import EnvInterpolation
import jaws_site.utils


conf = None


class ConfigurationError(Exception):
    pass


class Configuration(metaclass=jaws_site.utils.Singleton):

    """Configuration singleton class"""

    defaults = {
        "LOCAL_RPC_SERVER": {
            "host": "localhost",
            "port": "5672",
            "user": "guest",  # default from docker container
            "password": "guest",  # default from docker container
            "num_threads": 5,
            "max_retries": 5,
        },
        "CENTRAL_RPC_SERVER": {
            "host": "localhost",
            "port": "5672",
            "user": "guest",  # default from docker container
            "password": "guest",  # default from docker container
            "num_threads": 5,
            "max_retries": 5,
        },
        "CENTRAL_RPC_CLIENT": {
            "port": "5672",
        },
        "DB": {
            "host": "localhost",
            "port": "3306",
            "dialect": "mysql+mysqlconnector",
        },
    }
    required_params = {
        "LOCAL_RPC_SERVER": ["vhost"],
        "CENTRAL_RPC_SERVER": ["vhost"],
        "CENTRAL_RPC_CLIENT": ["host", "vhost", "user", "password"],
        "RUNS_ES_RPC_CLIENT": ["host", "vhost", "user", "password", "queue"],
        "PERFORMANCE_METRICS_ES_RPC_CLIENT": [
            "host",
            "vhost",
            "user",
            "password",
            "queue",
        ],
        "PERFORMANCE_METRICS": ["done_dir", "processed_dir"],
        "GLOBUS": [
            "client_id",
            "client_secret",
            "endpoint_id",
            "host_path",
        ],
        "DB": ["user", "password", "db"],
        "CROMWELL": ["url"],
        "SITE": ["id", "inputs_dir"],
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
        self.config = configparser.ConfigParser(interpolation=EnvInterpolation())
        self.config.read_dict(self.defaults)
        try:
            self.config.read(config_file)
        except configparser.ParsingError as error:
            logger.exception(f"Unable to load config file {config_file}: {error}")
            raise

        # validate config
        for section in self.required_params:
            if section not in self.config:
                error_msg = (
                    f"Config file, {config_file}, missing required section, {section}"
                )
                logger.error(error_msg)
                raise ValueError(error_msg)
            for key in self.required_params[section]:
                if key not in self.config[section]:
                    error_msg = f"Config file, {config_file}, missing required parameter, {section}/{key}"
                    logger.error(error_msg)
                    raise ValueError(error_msg)

        global conf
        conf = self

    def get(self, section: str, key: str, default=None) -> str:
        """Get a configuration value.

        :param section: name of config section
        :type section: str
        :param key: parameter key
        :type key: str
        :return: the value is always a string; typecast as necessary
        :rtype: str
        """
        if section not in self.config:
            raise ConfigurationError(f"Section {section} not defined in config obj")
        return self.config[section].get(key, default)

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
