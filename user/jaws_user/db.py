import logging
from urllib.parse import quote_plus
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import Column, String, Boolean
from sqlalchemy.ext.declarative import declarative_base
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

# ORM models

Base = declarative_base()


class User(Base):
    """Registered user"""

    __tablename__ = "users"
    id = Column(String(32), primary_key=True)
    email = Column(String(64), nullable=False)
    name = Column(String(64), nullable=True)
    is_admin = Column(Boolean, nullable=False, default=False)
    globus_id = Column(String(36), nullable=True)
    auth_refresh_token = Column(String(256), nullable=True)
    transfer_refresh_token = Column(String(256), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<User {self.id}>"


Base.metadata.create_all(engine)
