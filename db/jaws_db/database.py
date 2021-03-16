import logging
from urllib.parse import quote_plus
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, scoped_session


logger = logging.getLogger(__package__)


class DatabaseError(Exception):
    raise


class jaws_db:

    def __init__(self, dialect, user, password, host, port, db_name, pool_size=11, max_overflow=22, pool_recycle=3600, pool_pre_ping=True):
    self.db_name = db_name
    self.host = host
    self.pool_size = pool_size
    self.max_overflow = max_overflow
    self.pool_recycle = 3600
    self.pool_pre_ping = pool_pre_ping
    self.url = "%s://%s:%s@%s:%s/%s" % (
        "dialect"),
        "user"),
        quote_plus(password).encode('utf-8')),
        host,
        port,
        db_name)
    self.create_engine()

    def create_engine(self):
        logger.info(f"Connecting to db, {self.db_name} @ {self.host}")
        self.engine = None
        try:
            engine = create_engine(self.url, pool_size=self.pool_size, max_overflow=self.max_overflow, pool_recycle=self.pool_recycle, pool_timeout=self.pool_timeout, pool_pre_ping=self.pool_pre_ping)
        except Exception as error:
            err_msg = f"Unable to create engine for {self.db_name} @ {self.host}: {error}"
            logger.exception(err_msg)
            raise DatabaseError(err_msg)
        else:
            self.engine = engine

    def sessionmaker(self):
        if self.engine is None:
            self.create_engine()
        return sessionmaker(self.engine)

    def scoped_session(self):
        return scoped_session(self.sessionmaker())
