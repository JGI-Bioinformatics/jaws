import logging
import os
import configparser
from typing import Dict


MAX_SITE_ID_LEN = 8  # this matches the db column's varchar()


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


class ConfigItemNotFound(Exception):
    pass


class Configuration(metaclass=Singleton):
    """Configuration singleton class."""

    defaults = {
        "JAWS": {"name": "jaws", "version": "", "docs_url": ""},
        "DB": {"dialect": "mysql+mysqlconnector", "host": "localhost", "port": 3306},
        "RPC_SERVER": {
            "host": "localhost",
            "port": "5672",
            "queue": "central_rpc",
            "max_threads": 2,
            "max_retries": 3,
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
        "globus_endpoint",
        "input_dir",
        "output_dir",
        "max_ram_gb",
    ]
    default_site_params = {"globus_host_path": "/", "port": 5672, "message_ttl": 5}

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

        # validate config
        for section in self.required_params:
            if section not in self.config:
                error_msg = (
                    f"Config file, {config_file}, missing required section, {section}"
                )
                logger.error(error_msg)
                raise ConfigItemNotFound(error_msg)
            for key in self.required_params[section]:
                if key not in self.config[section]:
                    error_msg = f"Config file, {config_file}, missing required parameter, {section}/{key}"
                    logger.error(error_msg)
                    raise ConfigItemNotFound(error_msg)

        # init Sites
        self.sites = []
        for section in self.config.sections():
            if section.startswith("SITE:"):
                # extract site id
                site_id = section[len("SITE:") :].upper()
                if len(site_id) > MAX_SITE_ID_LEN:
                    # This matches the varchar() column type in the db
                    raise ConfigurationError(
                        f"Invalid Site ID: {site_id} (max. {MAX_SITE_ID_LEN} char.)"
                    )
                # save site id
                self.sites.append(site_id)
                # validate
                for key in self.required_site_params:
                    if key not in self.config[section]:
                        error_msg = f"Config file, {config_file}, missing required parameter, {section}/{key}"
                        logger.error(error_msg)
                        raise ContigItemNotFound(error_msg)

        # save singleton
        global conf
        conf = self

    def get(self, section: str, key: str) -> str:
        if section not in self.config:
            raise ConfigurationError(f"Section {section} not defined in config obj")
        if section in Configuration.defaults:
            return self.config[section].get(
                key, Configuration.defaults[section].get(key)
            )
        else:
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
        if site_id not in self.sites:
            raise ConfigItemNotFound
        section = f"SITE:{site_id}"
        return self.config[section].get(key, Configuration.default_site_params.get(key))

    def get_site_info(self, site_id: str) -> Dict[str, str]:
        """Returns public info about requested Site.

        :param site_id: The ID of the JAWS-Site
        :type site_id: str
        :return: Site parameters required to submit a run, if exists, None otherwise.
        :rtype: dict
        """
        site_id = site_id.upper()
        if site_id not in self.sites:
            raise ConfigItemNotFound
        section = f"SITE:{site_id}"
        # return only specific fields
        result = {"site_id": site_id}
        for key in [
            "globus_endpoint",
            "globus_host_path",
            "input_dir",
            "output_dir",
            "max_ram_gb",
        ]:
            default = Configuration.default_site_params.get(key, None)
            result[key] = self.config[section].get(key, default)
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
            raise ConfigItemNotFound
        section = f"SITE:{site_id}"

        # return only specific fields
        result = {"site_id": site_id}
        for key in [
            "host",
            "port",
            "user",
            "password",
            "vhost",
            "queue",
            "message_ttyl",
        ]:
            result[key] = get(
                self.config[section], key, Configuration.default_site_params.get(key)
            )
        return result

    def get_all_sites_rpc_params(self) -> Dict[str, Dict]:
        """Returns public info about requested Site.

        :return: Dict of site_id to dict of AMQP connection parameters.
        :rtype: dict
        """
        sites = {}
        for site_id in self.sites:
            sites[site_id] = self.get_site_rpc_params(site_id)
        return sites
