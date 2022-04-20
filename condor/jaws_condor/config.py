import logging
import os
import configparser


conf = None


class ConfigurationError(Exception):
    def __init__(self, message):
        super().__init__(message)


class Configuration():
    """Configuration singleton class."""

    config = None

    def __init__(self, config_file: str) -> None:
        """Constructor

        :param config_file: Path to config file in YAML format
        :type config_file: str
        """
        logger = logging.getLogger(__package__)
        logger.debug(f"Loading configuration from {config_file}")
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"{config_file} does not exist")
        self.config = configparser.ConfigParser()
        try:
            self.config.read(config_file)
        except Exception as error:
            logger.exception(f"Unable to load config file {config_file}: {error}")
            raise

        # save singleton
        global conf
        conf = self

    def get_site_id(self):
        """Returns multiprocessing connection info.

        :return: Multiprocessing parameters.
        :rtype: dict
        """
        section = "SITE"
        s = self.config[section]
        params = {
            "site_id":  s["site_id"],
        }
        return s["site_id"]