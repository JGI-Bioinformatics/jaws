import logging
from urllib.parse import quote_plus
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from jaws_central import models


class Singleton(type):

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class JawsDb(metaclass=Singleton):
    """Database singleton class"""
    engine = None
    session = None

    def __init__(self, conf):
        """Constructor

        :param conf: An initialized JawsConfig object
        :type conf: obj
        :return: JawsDb singleton object
        :rtype: obj
        """
        logger = logging.getLogger(__package__)
        logger.info("Initializing db connection")
        url = "%s://%s:%s@%s/%s" % (
            conf.get("db", "dialect"),
            conf.get("db", "user"),
            quote_plus(conf.get("db", "password")),
            conf.get("db", "host"),
            conf.get("db", "db"))
        self.engine = create_engine(url, pool_size=3, pool_recycle=3600, pool_pre_ping=True)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        models.create_all(self.engine)

    def engine(self):
        if self.engine is None:
            raise Exception("Db not initialized; run init_db first")
        return self.engine

    def session(self):
        if self.engine is None:
            raise Exception("Db not initialized; run init_db first")
        return self.session


db = None
