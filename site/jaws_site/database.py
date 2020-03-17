import logging
from urllib.parse import quote_plus
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from jaws_site import models


db = None


class Singleton(type):

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class DatabaseError(Exception):
    pass


class JawsDb(metaclass=Singleton):
    """Database singleton class"""
    engine = None

    def __init__(self, conf):
        """Constructor

        :param conf: An initialized JawsConfig object
        :type conf: obj
        :return: JawsDb singleton object
        :rtype: obj
        """
        logger = logging.getLogger(__package__)
        logger.debug("Initializing db connection")
        url = "%s://%s:%s@%s/%s" % (
            conf.get("DB", "dialect"),
            conf.get("DB", "user"),
            quote_plus(conf.get("DB", "password").encode('utf-8')),
            conf.get("DB", "host"),
            conf.get("DB", "db"))
        try:
            self.engine = create_engine(url, pool_size=3, pool_recycle=3600, pool_pre_ping=True)
        except Exception as e:
            raise DatabaseError(f"Unable to connect to db: {e}")
        models.Base.metadata.create_all(self.engine)
        global db
        db = self

    def engine(self):
        if self.engine is None:
            raise DatabaseError("Db not initialized; run init_db first")
        return self.engine

    def session(self):
        """Return a new session obj."""
        if self.engine is None:
            raise DatabaseError("Db not initialized; run init_db first")
        Session = sessionmaker(bind=self.engine)
        session = Session()
        return session
