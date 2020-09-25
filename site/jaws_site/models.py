"""
SQLAlchemy models for persistent data structures.
"""

from datetime import datetime
from sqlalchemy import Column, DateTime, String, Integer, Boolean, ForeignKey
from jaws_site.database import Base


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
    reason = Column(String(1024), nullable=True)
    sent = Column(Boolean, default=False, nullable=False)


def create_all(engine, session):
    """Create all tables if not exist"""
    Base.metadata.create_all(engine)
    session.commit()
