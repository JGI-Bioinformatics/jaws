import logging
from urllib.parse import quote_plus
from sqlalchemy import create_engine, MetaData
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import TimeoutError
from jaws_central import config


logger = logging.getLogger(__package__)

conf = config.Configuration()
params = conf.get_section("DB")
logger.info(f"Connecting to db, {params.get('db')} @ {params.get('host')}")
url = "%s://%s:%s@%s:%s/%s" % (
    params.get("dialect"),
    params.get("user"),
    quote_plus(params.get("password").encode('utf-8')),
    params.get("host"),
    params.get("port"),
    params.get("db"))
try:
    engine = create_engine(url, pool_size=5, max_overflow=10, pool_recycle=3600, pool_timeout=30, pool_pre_ping=True)
except TimeoutError as error:
    logger.exception(f"Db unavailable: {error}")
    raise
except Exception as error:
    logger.exception(error)
    raise

Session = sessionmaker(bind=engine)

metadata = MetaData()
Base = declarative_base(metadata=metadata)
