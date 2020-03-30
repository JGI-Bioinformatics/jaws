import logging
import os
import configparser
from typing import Dict


conf = None


class Singleton(type):

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class ConfigurationError(Exception):
    def __init__(self, message):
        super().__init__(message)


class JawsConfig(metaclass=Singleton):
    """Configuration singleton class.
    """

    config = None

    def __init__(self, config_file=None):
        """Constructor

        :param config_file: Path to config file in YAML format
        :type config_file: str
        """
        logger = logging.getLogger(__package__)
        logger.debug(f"Loading configuration from {config_file}")
        if not config_file:
            raise FileNotFoundError("config file not specified")
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"{config_file} does not exist")
        self.config = configparser.ConfigParser()
        self.config.read(config_file)
        required_sections = ("DB", "GLOBUS")
        for section in required_sections:
            if section not in self.config:
                error_msg = f'Config file, {config_file}, missing required "{section}" section'
                logger.warn(error_msg)
                raise ConfigurationError(error_msg)
        global conf
        conf = self
        self._init_sites()

    def _init_sites(self):
        self.sites = {}
        for section in self.config.sections():
            if section.startswith("SITE:"):
                site_id = section[len("SITE:"):].upper()
                if len(site_id) > 8:
                    raise ConfigurationError(f"Invalid Site ID: {site_id} (max. 8 char.)")
                s = self.sites[site_id] = self.config[section]
                s["staging_dir"] = os.path.join(s["globus_basepath"], s["staging_subdir"])

    def get(self, section: str, key: str) -> str:
        if section not in self.config:
            raise ConfigurationError(f"Section {section} not defined in config obj")
        return self.config[section].get(key)

    def get_site(self, site_id: str, key: str) -> str:
        """Retrieve Site config parameter; syntactic sugar.

        :param site_id: Unique ID of a JAWS-Site
        :type site_id: str
        :param key: The desired parameter
        :type key: str
        :return: The configuration value
        :rtype: str
        """
        site_id = site_id.upper()
        section = f"SITE:{site_id}"
        return self.get(section, key)

    def get_site_info(self, site_id: str) -> Dict[str, str]:
        """Returns public info about requested Site.

        :param site_id: The ID of the JAWS-Site
        :type site_id: str
        :return: Site parameters required to submit a run, if exists, None otherwise.
        :rtype: dict
        """
        site_id = site_id.upper()
        if site_id not in self.sites:
            return None
        section = f"SITE:{site_id}"
        s = self.config[section]
        result = {
            "site_id": site_id,
            "globus_endpoint": s["globus_endpoint"],
            "globus_basepath": s["globus_basepath"],
            "staging_subdir": s["staging_subdir"],
            "max_ram_gb": s["max_ram_gb"],
        }
        return result
