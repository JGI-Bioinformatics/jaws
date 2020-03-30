import os
import configparser
import jaws_site.utils


USER_JAWSRC = os.path.join(os.path.expanduser("~"), ".jawsrc", "jaws-site.ini")
JAWS_SITE_CONFIG = "JAWS_SITE_CONFIG"


class Configuration(metaclass=jaws_site.utils.Singleton):

    """Configuration singleton class"""
    defaults = {"AMQP": {
        "host": "localhost",
        "vhost": "/",
        "user": "guest",  # default from docker container
        "password": "guest",  # default from docker container
        "queue": "test"
        },
        "RPC": {
            "host": "localhost",
            "vhost": "/",
            "user": "jaws",  # set in docker-compose
            "password": "jawstest",  # set in docker-compose
            "queue": "test",
            "num_threads": 5,
            "max_retries": 5,
        },
        "GLOBUS": {
            "client_id": "",
            "endpoint_id": "jaws-testing",
            "root_dir": "",
            "subdir": ""
        },
        "DB": {
            "host": "localhost",
            "port": "3306",
            "user": "jaws",
            "password": "jawstest",
            "db": "jaws-workflows",
            "dialect": "mysql+mysqlconnector"
        },
        "CROMWELL": {
            "host": "localhost",
            "port": "8000",
            "workflows_url": "/api/workflows/v2",
            "engine_status_url": "/engine/v2"
        },
        "SITE": {
            "id": "local",
            "staging_subdirectory": "$HOME/staging",
            "results_subdirectory": "$HOME/results"
        }
    }

    def __init__(self, config_file: str = None):
        """Constructor sets global singleton.

        :param config_file: Path to configuration file in INI format
        :type config_file: str
        """
        self.config = configparser.ConfigParser()
        self.config.read_dict(self.defaults)

        if config_file and not os.path.exists(config_file):
            raise FileNotFoundError(f"{config_file} does not exist")
        else:
            self.config.read(config_file)

    def get(self, section: str, key: str) -> str:
        """Get a configuration value.

        :param section: name of config section
        :type section: str
        :param key: parameter key
        :type key: str
        :return: the value is always a string; typecast as necessary
        :rtype: str
        """
        return self.config.get(section, key)


def get_config():
    """Find config file by checking env var, cwd, and home, in that order.

    :param config_file: filename of configuration file
    :type config_file: str
    :param env_var: name of config file environment variable
    :type env_var: str
    :return: absolute path if found, None otherwise.
    :rtype: str
    """
    path = ""
    if JAWS_SITE_CONFIG in os.environ:
        path = os.environ.get(JAWS_SITE_CONFIG)
    elif os.path.exists(USER_JAWSRC):
        path = USER_JAWSRC
    return Configuration(path)


jaws_config = get_config()
