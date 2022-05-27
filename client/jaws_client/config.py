"""
Configuration singleton, loads values from provided INI files.
"""

import logging
import os
import configparser
import shutil


JAWS_CONFIG_ENV = "JAWS_CLIENT_CONFIG"
JAWS_USER_CONFIG_ENV = "JAWS_USER_CONFIG"


class ConfigurationError(Exception):
    pass


class Configuration:
    """Configuration class variables should be initialized at the start of the program."""

    defaults = {
        "USER": {"token": ""},
        "JAWS": {
            "name": "JAWS",
            "site_id": "",
            "url": "http://localhost:5000",
            "womtool_jar": "",
            "inputs_dir": "/tmp",
            "downloads_dir": "/tmp",
        },
    }

    required_jaws_params = {
        "JAWS": [
            "name",
            "site_id",
            "url",
            "womtool_jar",
            "inputs_dir",
            "downloads_dir",
        ],
    }
    required_user_params = {"USER": ["token"]}

    config = None

    def __init__(
        self, jaws_config_file: str = None, user_config_file: str = None
    ) -> None:
        """Initialize the configuration object.  Use class variables if paths not provided.

        :param jaws_config_file: Path to configuration file of JAWS deployment
        :type jaws_config_file: str
        :param user_config_file: Path to configuration file of user info
        :type user_config_file: str
        """
        if self.config is not None:
            return
        if not jaws_config_file:
            jaws_config_file = os.environ.get(JAWS_CONFIG_ENV)
            if not jaws_config_file:
                raise ValueError("config path required")
        if not user_config_file:
            user_config_file = os.environ.get(JAWS_USER_CONFIG_ENV)
            if not user_config_file:
                raise ValueError("user config path required")
        logger = logging.getLogger(__package__)
        logger.debug(f"Loading config from {jaws_config_file}, {user_config_file}")
        if not os.path.isfile(jaws_config_file):
            raise FileNotFoundError(f"{jaws_config_file} does not exist")
        if not os.path.isfile(user_config_file):
            raise FileNotFoundError(f"{user_config_file} does not exist")
        self.config = configparser.ConfigParser()
        self.config.read_dict(self.defaults)
        try:
            self.config.read(jaws_config_file)
        except configparser.Error as error:
            logger.exception(f"Invalid config file, {jaws_config_file}: {error}")
            raise

        # validate config
        for section in self.required_jaws_params:
            if section not in self.config:
                error_msg = f"JAWS config missing required section, {section}"
                logger.error(error_msg)
                raise ValueError(error_msg)
            for key in self.required_jaws_params[section]:
                if key not in self.config[section]:
                    error_msg = (
                        f"JAWS config missing required parameter, {section}/{key}"
                    )
                    logger.error(error_msg)
                    raise ValueError(error_msg)

        # validate womtool jar path
        womtool_path = self.config["JAWS"]["womtool_jar"]
        if womtool_path:
            if not os.path.isfile(womtool_path):
                raise FileNotFoundError(f"womtool jar does not exist: {womtool_path}")
            self.config["JAWS"][
                "womtool"
            ] = f"java -jar {self.config['JAWS']['womtool_jar']}"
        else:
            womtool_path = shutil.which("womtool")
            if not womtool_path:
                raise FileNotFoundError("womtool not found in path or config file")
            self.config["JAWS"]["womtool"] = womtool_path

        # load user config and copy into config obj
        user_config = configparser.ConfigParser()
        try:
            user_config.read(user_config_file)
        except configparser.Error as error:
            logger.exception(f"Invalid config file, {user_config_file}: {error}")
            raise ConfigurationError(f"Invalid config file, {user_config_file}: {error}")
        for section in self.required_user_params:
            self.config[section] = {}
            if section not in user_config:
                error_msg = f"JAWS config missing required section, {section}"
                logger.error(error_msg)
                raise ValueError(error_msg)
            for key in self.required_user_params[section]:
                if key not in user_config[section]:
                    error_msg = (
                        f"JAWS config missing required parameter, {section}/{key}"
                    )
                    logger.error(error_msg)
                    raise ValueError(error_msg)
                self.config[section][key] = user_config[section][key]

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
