import logging
from datetime import datetime
from urllib.parse import quote_plus
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy import Column, DateTime, String, Integer
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
        url, pool_size=11, max_overflow=22, pool_recycle=3600, pool_pre_ping=True
    )
except Exception as error:
    logger.exception(error)
    raise

Session = scoped_session(sessionmaker(bind=engine))
Session.configure(bind=engine)

Base = declarative_base()


class Job_Log(Base):
    """
    Log state transitions log.
    Initially, records are inserted when a state transition log is received from JTM;
    however JTM doesn't know the run_id, task_name, or attempt so they are NULL.
    The Daemon process will fill in those fields by querying the db and cromwell;
    until then, they will not be found by selecting run_id.
    """

    __tablename__ = "job_logs"
    run_id = Column(Integer, nullable=True)
    cromwell_run_id = Column(String(36), nullable=False)
    cromwell_job_id = Column(Integer, primary_key=True)
    task_name = Column(String(128), nullable=True)
    attempt = Column(Integer, nullable=True)
    status_from = Column(String(32), primary_key=True)
    status_to = Column(String(32), primary_key=True)
    timestamp = Column(DateTime, nullable=False, default=datetime.utcnow)
    reason = Column(String(1024), nullable=True)


def create_all(engine, session):
    """Create all tables if not exist"""
    Base.metadata.create_all(engine)
    session.commit()
