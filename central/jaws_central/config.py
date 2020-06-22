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

    def _destructor(cls):
        del cls._instances[cls]


class ConfigurationError(Exception):
    def __init__(self, message):
        super().__init__(message)


class Configuration(metaclass=Singleton):
    """Configuration singleton class."""

    defaults = {
        "DB": {
            "dialect": "mysql+mysqlconnector",
            "host": "localhost",
            "port": 3306,
            "user": "jaws",
            "password": "",
            "db": "jaws",
        },
        "GLOBUS": {"client_id": ""},
        "RPC_SERVER": {
            "host": "localhost",
            "port": "5672",
            "user": "jaws_rmq_user",
            "password": "jaws_rmq_pass",
            "vhost": "jaws_central",
            "max_threads": 5,
            "max_retries": 3
        }
    }

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
        self.config.read_dict(self.defaults)
        try:
            self.config.read(config_file)
        except Exception as error:
            logger.exception(f"Unable to load config file {config_file}: {error}")
            raise
        global conf
        conf = self
        self._init_sites()

    def _init_sites(self):
        self.sites = {}
        for section in self.config.sections():
            if section.startswith("SITE:"):
                site_id = section[len("SITE:"):].upper()
                if len(site_id) > 8:
                    raise ConfigurationError(
                        f"Invalid Site ID: {site_id} (max. 8 char.)"
                    )
                s = self.sites[site_id] = self.config[section]
                s["staging_dir"] = os.path.join(
                    s["globus_basepath"], s["staging_subdir"]
                )

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

    def get_site_rpc_params(self, site_id: str) -> Dict[str, str]:
        """Returns AMQP connection info for the site.

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
        params = {
            "host": s["host"],
            "port": s["port"],
            "user": s["user"],
            "password": s["password"],
            "vhost": s["vhost"]
        }
        return params

    def get_all_sites_rpc_params(self) -> Dict[str, Dict]:
        """Returns public info about requested Site.

        :return: Dict of site_id to dict of AMQP connection parameters.
        :rtype: dict
        """
        sites = {}
        for site_id in self.sites:
            sites[site_id] = self.get_site_rpc_params(site_id)
        return sites
