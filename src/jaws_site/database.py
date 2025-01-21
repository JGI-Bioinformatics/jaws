import logging
from urllib.parse import quote_plus

import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import registry, sessionmaker

from sqlalchemy.schema import Index
from sqlalchemy import inspect
from sqlalchemy.exc import SQLAlchemyError
from jaws_site import config

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
        pool_size=11,
        max_overflow=22,
        pool_recycle=3600,
        pool_timeout=30,
        pool_pre_ping=True,
    )
except Exception as error:
    logger.exception(error)
    raise

Session = sessionmaker(engine, autobegin=True, autocommit=False)

mapper_registry = registry()
Base = mapper_registry.generate_base()


def create_index_if_not_exists(index_name: str, table_name: str, col_name: str, engine: sqlalchemy.engine.base.Engine):

    """Create an index if it does not already exist.

    :param index_name: name of the index
    :type index_name: str
    :param table_name: name of the table
    :type table_name: str
    :param col_name: name of the column
    :type col_name: str
    :param engine: SQLAlchemy engine
    :type engine: sqlalchemy.engine.base.Engine
    """

    try:
        inspector = inspect(engine)
        indexes = inspector.get_indexes(table_name)
        if index_name not in [index["name"] for index in indexes]:
            logger.debug(f"Creating index {index_name} on {table_name}.{col_name}")
            Index(index_name, getattr(Base.metadata.tables[table_name].c, col_name)).create(engine)
        else:
            logger.debug(f"Index {index_name} already exists on {table_name}.{col_name}")
    except SQLAlchemyError as e:
        logger.error(f"Failed to create index {index_name} on {table_name}.{col_name}: {e}")
