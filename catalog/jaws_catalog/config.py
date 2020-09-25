import logging
import os
import configparser


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
        "CATALOG": {"port": "6000"}
        "DB": {"dialect": "mysql+mysqlconnector", "host": "localhost", "port": 3306},
    }
    required_params = {
        "DB": ["user", "password", "db"],
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

        # save singleton
        global conf
        conf = self

    def get(self, section: str, key: str) -> str:
        if section not in self.config:
            raise ConfigurationError(f"Section {section} not defined in config obj")
        return self.config[section].get(key)

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
