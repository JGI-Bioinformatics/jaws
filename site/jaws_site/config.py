"""
Configuration singleton, loads values from provided INI infile.
"""

import logging
import os
import configparser
from urllib.parse import quote_plus
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from jaws_site import models


class Singleton(type):

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class JawsConfig(metaclass=Singleton):
    config = None
    db = None
    session = None

    def __init__(self, config_file=None):
        """Constructor

        :param config_file: Path to configuration file in INI format
        :type config_file: str
        """
        logger = logging.getLogger(__package__)
        logger.debug('loading configuration...')
        if not config_file:
            raise FileNotFoundError("config file not specified")
        if not os.path.isfile(config_file):
            raise FileNotFoundError(f"{config_file} does not exist")
        self.config = configparser.ConfigParser()
        self.config.read(config_file)

    def get_db(self, value):
        return self.config.get("DB", value)

    def get_amqp(self, value):
        return self.config.get("AMQP", value)

    def get_rpc(self, value):
        return self.config.get("RPC", value)

    def get_globus(self, value):
        return self.config.get("GLOBUS", value)

    def get_cromwell(self, value):
        return self.config.get("CROMWELL", value)

    def get_site(self, value):
        return self.config.get("SITE", value)

    def init_db(self):
        """SQLAlchemy connection to RDB
        """
        url = "%s://%s:%s@%s/%s" % (
            self.get_db("dialect"),
            self.get_db("user"),
            quote_plus(self.get_db("password")),
            self.get_db("host"),
            self.get_db("db"))
        self.db = create_engine(url)
        Session = sessionmaker(bind=self.db)
        self.session = Session()
        models.create_all()

    def db(self):
        if self.db is None:
            raise Exception("db not initialized; run init_db() first")
        return self.db

    def session(self):
        if self.db is None:
            raise Exception("db not initialized; run init_db() first")
        return self.session
