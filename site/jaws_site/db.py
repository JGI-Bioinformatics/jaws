import logging
from datetime import datetime
from urllib.parse import quote_plus
from sqlalchemy import create_engine, Column, DateTime, String, Integer, Boolean, ForeignKey, UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, scoped_session
from jaws_site import config

logger = logging.getLogger(__package__)

params = config.conf.get_section("DB")
logger.info(f"Connecting to db, {params.get('db')} @ {params.get('host')}")
url = "%s://%s:%s@%s:%s/%s" % (
    params.get("dialect"),
    params.get("user"),
    quote_plus(params.get("password").encode('utf-8')),
    params.get("host"),
    params.get("port"),
    params.get("db"))
try:
    engine = create_engine(url, pool_size=15, max_overflow=30, pool_recycle=3600, pool_timeout=30, pool_pre_ping=True)
except Exception as error:
    logger.exception(error)
    raise

Session = scoped_session(sessionmaker(bind=engine))
Session.configure(bind=engine)

Base = declarative_base()


class Run(Base):
    """
    Analysis runs are the execution of workflows on specific inputs.
    Only the active runs are stored; finished/failed runs are deleted.
    Run info is permanently recorded in Central's "runs" table instead.
    """

    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    submission_id = Column(String(36), nullable=False)
    cromwell_run_id = Column(String(36), nullable=True)
    cromwell_workflow_dir = Column(String(4096), nullable=True)
    status = Column(String(32), nullable=False)
    user_id = Column(String(32), nullable=False)
    submitted = Column(DateTime, nullable=False, default=datetime.utcnow)
    updated = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    upload_task_id = Column(String(36), nullable=False)
    output_endpoint = Column(String(36), nullable=False)
    output_dir = Column(String(256), nullable=False)
    download_task_id = Column(String(36), nullable=True)
    email = Column(String(64), nullable=False)
    transfer_refresh_token = Column(String(256), nullable=False)


class Run_Log(Base):
    """
    Run state transitions log.
    """

    __tablename__ = "run_logs"
    run_id = Column(Integer, ForeignKey("runs.id"), primary_key=True)
    status_from = Column(String(32), primary_key=True)
    status_to = Column(String(32), primary_key=True)
    timestamp = Column(DateTime, nullable=False, default=datetime.utcnow)
    reason = Column(String(1024), nullable=False, default="")
    sent = Column(Boolean, default=False, nullable=False)


class Job_Log(Base):
    """
    Log state transitions log.
    Initially, records are inserted when a state transition log is received from JTM;
    however JTM doesn't know the run_id, task_name, or attempt so they are NULL.
    The Daemon process will fill in those fields by querying the db and cromwell;
    until then, they will not be found by selecting run_id.
    """

    __tablename__ = "job_logs"
    run_id = Column(Integer, ForeignKey("runs.id"), nullable=True)
    cromwell_run_id = Column(String(36), nullable=False)
    cromwell_job_id = Column(Integer, primary_key=True)
    task_name = Column(String(128), nullable=True)
    attempt = Column(Integer, nullable=True)
    status_from = Column(String(32), primary_key=True)
    status_to = Column(String(32), primary_key=True)
    timestamp = Column(DateTime, nullable=False, default=datetime.utcnow)
    reason = Column(String(1024), nullable=False, default="")


class Task_Log(Base):
    """
    Task state transitions log.
    """

    __tablename__ = "task_logs"
    id = Column(Integer, primary_key=True)
    cromwell_run_id = Column(String(36), nullable=False)
    task_id = Column(Integer, nullable=False)
    status_from = Column(String(32), nullable=False)
    status_to = Column(String(32), nullable=False)
    timestamp = Column(DateTime, nullable=False, default=datetime.utcnow)
    reason = Column(String(1024), nullable=False, default="")
    UniqueConstraint("task_id", "status_from")


def create_all(engine, session):
    """Create all tables if not exist"""
    Base.metadata.create_all(engine)
    session.commit()
