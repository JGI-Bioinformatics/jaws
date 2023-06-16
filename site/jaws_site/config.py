import os
import logging
import configparser
from jaws_site.env_interpolation import EnvInterpolation, JAWSConfigParser
import jaws_site.utils


conf = None
DEFAULT_MAX_RAM_GB = 56
DEFAULT_MAX_CPU = 32


class ConfigurationError(Exception):
    pass


class Configuration(metaclass=jaws_site.utils.Singleton):

    """Configuration singleton class"""

    defaults = {
        "DB": {
            "host": "localhost",
            "port": "3306",
            "dialect": "mysql+mysqlconnector",
        },
        "RMQ": {
            "host": "localhost",
            "port": "5672",
            "user": "guest",  # default from docker container
            "password": "guest",  # default from docker container
            "num_threads": 5,
            "max_retries": 3,
        },
        "SITE": {
            "max_user_active_runs": 0,
            "max_transfer_threads": 7,
            "file_permissions": 755,
        },
    }
    required_params = {
        "PERFORMANCE_METRICS": [
            "done_dir",
            "processed_dir",
            "running_dir",
            "cleanup_time",
        ],
        "GLOBUS": [
            "endpoint_id",
            "host_path",
        ],
        "DB": ["user", "password"],
        "CROMWELL": ["url"],
        "SITE": ["id", "deployment", "inputs_dir"],
    }

    config = None

    def __init__(self, config_file: str, env_prefix: str = None) -> None:
        """Constructor sets global singleton.

        :param config_file: Path to configuration file in INI format
        :type config_file: str
        """
        logger = logging.getLogger(__package__)
        logger.debug(f"Loading config from {config_file}")
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"{config_file} does not exist")
        self.config = JAWSConfigParser(
            interpolation=EnvInterpolation(), env_override=env_prefix
        )
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
        if section not in self.config:
            error_msg = f"Config missing requested section: {section}"
            raise ValueError(error_msg)
        sect_conf = self.config[section]
        for key, value in sect_conf.items():
            result[key] = value
        return result

    def get_site_config(self):
        result = {}
        result["max_ram_gb"] = self.get("SITE", "max_ram_gb", DEFAULT_MAX_RAM_GB)
        result["max_cpu"] = self.get("SITE", "max_cpu", DEFAULT_MAX_CPU)
        result["inputs_dir"] = self.get("SITE", "inputs_dir")
        result["downloads_dir"] = self.get("SITE", "downloads_dir")
        result["access_group"] = self.get("SITE", "access_group")
        result["globus_host_path"] = self.get("GLOBUS", "host_path")
        result["globus_endpoint"] = self.get("GLOBUS", "endpoint_id")
        result["file_sync_delay_sec"] = self.get("SITE", "file_sync_delay_sec", 0)
        return result
