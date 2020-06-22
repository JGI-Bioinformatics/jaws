"""
SQLAlchemy models for persistent data structures.
"""

from datetime import datetime
from sqlalchemy import (
    Column,
    DateTime,
    String,
    Integer,
    ForeignKey
)
from jaws_site.database import Base


class Run(Base):
    """Analysis runs are the execution of workflows on specific inputs."""

    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    submission_id = Column(String(36), nullable=False)
    cromwell_run_id = Column(String(36), nullable=True)
    status = Column(String(16), nullable=False)
    user_id = Column(String(16), nullable=False)
    site_id = Column(String(8), nullable=False)
    submitted = Column(DateTime, nullable=False, default=datetime.utcnow)
    updated = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    upload_task_id = Column(String(36), nullable=False)
    output_endpoint = Column(String(36), nullable=False)
    output_dir = Column(String(256), nullable=False)
    download_task_id = Column(String(36), nullable=True)
    email = Column(String(64), nullable=False)
    transfer_refresh_token = Column(String(256), nullable=False)


class Run_Log(Base):
    """Run state transitions log"""

    __tablename__ = "run_logs"
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey("runs.id"), nullable=False)
    status_from = Column(String(32), nullable=False)
    status_to = Column(String(32), nullable=False)
    timestamp = Column(DateTime, nullable=False, default=datetime.utcnow)
    reason = Column(String(1024), nullable=True)


class Job_Log(Base):
    """Job state transitions are recorded only until sent to Central"""

    __tablename__ = "job_logs"
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey("runs.id"), nullable=False)
    cromwell_job_id = Column(Integer, nullable=True)
    task_name = Column(String(128), nullable=True)
    attempt = Column(Integer, nullable=True)
    status_from = Column(String(32), nullable=False)
    status_to = Column(String(32), nullable=False)
    timestamp = Column(DateTime, nullable=False, default=datetime.utcnow)
    reason = Column(String(1024), nullable=True)


def create_all(engine):
    """Create all tables if not exist"""
    Base.metadata.create_all(engine)
