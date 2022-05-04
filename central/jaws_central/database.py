import logging
from urllib.parse import quote_plus
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, scoped_session
from jaws_central import config


logger = logging.getLogger(__package__)

params = config.conf.get_section("DB")
logger.info(f"Connecting to db, {params.get('db')} @ {params.get('host')}")
url = "%s://%s:%s@%s:%s/%s" % (
    params.get("dialect"),
    params.get("user"),
    quote_plus(params.get("password").encode("utf-8")),
    params.get("host"),
    params.get("port"),
    params.get("db"),
)
try:
    engine = create_engine(
        url,
        pool_size=5,
        max_overflow=10,
        pool_recycle=3600,
        pool_timeout=30,
        pool_pre_ping=True,
    )
except Exception as error:
    logger.exception(error)
    raise

Session = scoped_session(sessionmaker(bind=engine))
Session.configure(bind=engine)

Base = declarative_base()
