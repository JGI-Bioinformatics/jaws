"""
SQLAlchemy models for persistent data structures.
"""

from datetime import datetime
from sqlalchemy import Column, DateTime, String, Integer
from jaws_site.database import Base


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
