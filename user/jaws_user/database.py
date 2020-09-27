import logging
from urllib.parse import quote_plus
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from jaws_user import config


logger = logging.getLogger(__package__)

# get db config
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

# init db engine
try:
    engine = create_engine(url, pool_size=3, pool_recycle=3600, pool_pre_ping=True)
except Exception as error:
    logger.exception(error)
    raise

# init sessionmaker
Session = sessionmaker(bind=engine)
