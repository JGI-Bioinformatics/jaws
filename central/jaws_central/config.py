import logging
import os
import configparser
from jaws_central.env_interpolation import EnvInterpolation
from typing import Dict


DEFAULT_AMQP_PORT = 5672
DEFAULT_RPC_MESSAGE_TTL = 5
MAX_SITE_ID_LEN = 8  # this matches the MYSQL column's varchar()


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
        "JAWS": {"name": "jaws", "version": "", "docs_url": ""},
        "DB": {"dialect": "mysql+mysqlconnector", "host": "localhost", "port": 3306},
        "RPC_SERVER": {
            "host": "localhost",
            "port": "5672",
            "user": "guest",  # default from docker container
            "password": "guest",  # default from docker container
            "queue": "central_rpc",
            "num_threads": 5,
            "max_retries": 5,
        },
        "HTTP": {"auth_url": "localhost", "auth_port": "3000", "rest_port": "5000"},
    }
    required_params = {
        "DB": ["user", "password", "db"],
        "GLOBUS": ["client_id", "client_secret"],
        "RPC_SERVER": ["user", "password", "vhost"],
    }
    required_site_params = [
        "host",
        "user",
        "password",
        "vhost",
        "globus_host_path",
        "globus_endpoint",
        "inputs_dir",
        "max_ram_gb",
    ]

    config = None

    def __init__(self, config_file: str, env_prefix: str = None) -> None:
        """Constructor

        :param config_file: Path to config file in YAML format
        :type config_file: str
        """
        logger = logging.getLogger(__package__)
        logger.debug(f"Loading configuration from {config_file}")
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"{config_file} does not exist")
        self.config = configparser.ConfigParser(interpolation=EnvInterpolation(env_override=env_prefix))
        self.config.read_dict(self.defaults)
        try:
            self.config.read(config_file)
        except Exception as error:
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

        # init Sites
        self.sites = {}
        for section in self.config.sections():
            if section.startswith("SITE:"):
                site_id = section[len("SITE:") :].upper()  # noqa
                if len(site_id) > MAX_SITE_ID_LEN:
                    raise ConfigurationError(
                        f"Invalid Site ID: {site_id} (max. {MAX_SITE_ID_LEN} char.)"
                    )
                # validate
                for key in self.required_site_params:
                    if key not in self.config[section]:
                        error_msg = f"Config file, {config_file}, missing required parameter, {section}/{key}"
                        logger.error(error_msg)
                        raise ValueError(error_msg)
                # copy
                self.sites[site_id] = {}
                for key in self.config[section]:
                    self.sites[site_id][key] = self.config[section][key]
                self.sites[site_id]["inputs_dir"] = os.path.join(
                    self.sites[site_id]["globus_host_path"],
                    self.sites[site_id]["inputs_dir"],
                )

        # save singleton
        global conf
        conf = self

    def get(self, section: str, key: str, default=None) -> str:
        if section not in self.config:
            raise ConfigurationError(f"Section {section} not defined in config obj")
        return self.config[section].get(key, default)

    def get_site(self, site_id: str) -> dict:
        """Retrieve Site config section
        :param site_id: Unique ID of a JAWS-Site
        :type site_id: str
        :return: All of a Site's configuration parameters
        :rtype: dict
        """
        return self.get_section(f"SITE:{site_id}")

    def get_site_param(self, site_id: str, key: str) -> str:
        """Retrieve Site config parameter; syntactic sugar.

        :param site_id: Unique ID of a JAWS-Site
        :type site_id: str
        :param key: The desired parameter
        :type key: str
        :return: The configuration value
        :rtype: str
        """
        site_id = site_id.upper()
        return self.sites[site_id].get(key)

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
            "globus_host_path": s["globus_host_path"],
            "inputs_dir": s["inputs_dir"],
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
            "port": int(s.get("port", DEFAULT_AMQP_PORT)),
            "user": s["user"],
            "password": s["password"],
            "vhost": s["vhost"],
            "queue": s["queue"],
            "message_ttl": int(s.get("message_ttl", DEFAULT_RPC_MESSAGE_TTL)),
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
