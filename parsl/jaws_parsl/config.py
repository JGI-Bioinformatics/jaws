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

    def get_rmq_params(self):
        """Returns RMQ connection info.

        :return: RabbitMQ parameters.
        :rtype: dict
        """
        section = f"RMQ"
        s = self.config[section]
        params = {
            "user": s["user"],
            "password": s["password"],
            "host": s["host"],
            "vhost": s["vhost"],
            "port": int(s["port"]),
        }
        print(params)
        return params

    def get_rpc_params(self):
        """Returns site RPC connection info.

        :return: RPC parameters.
        :rtype: dict
        """
        section = f"SITE_RPC_CLIENT"
        s = self.config[section]
        params = {
            "user": s["user"],
            "password": s["password"],
            "host": s["host"],
            "vhost": s["vhost"],
            "port": int(s["port"]),
            "queue": s["queue"]
        }
        print(params)
        return params

