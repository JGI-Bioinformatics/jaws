"""
SQLAlchemy models for persistent data structures.
"""

from datetime import datetime
from sqlalchemy import Column, DateTime, String, Integer, Boolean, ForeignKey
from sqlalchemy.dialects.mysql import MEDIUMTEXT, LONGTEXT
from jaws_site.database import Base


def same_as(column_name: str):
    """Function sets the default value of a column to the value in another column.

    :param column_name: name of the column
    :type column_name: str
    :return: function which retrieves the value of the specified column
    :rtype: function
    """

    def default_function(context):
        return context.current_parameters.get(column_name)

    return default_function


class Run(Base):
    """
    Analysis runs are the execution of workflows on specific inputs.
    Only the active runs are stored; finished/failed runs are deleted.
    Run info is permanently recorded in Central's "runs" table instead.
    """

    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    submission_id = Column(String(36), nullable=False)
    input_site_id = Column(String(16), nullable=False)
    caching = Column(Boolean, nullable=False, default=True)
    cromwell_run_id = Column(String(36), nullable=True)
    status = Column(String(32), nullable=False)
    result = Column(String(9), nullable=True)
    user_id = Column(String(32), nullable=False)
    submitted = Column(DateTime, nullable=False, default=datetime.utcnow)
    updated = Column(
        DateTime,
        default=same_as("submitted"),
        onupdate=datetime.utcnow,
    )
    workflow_root = Column(String(1024), nullable=True)
    workflow_name = Column(String(64), nullable=True)


class Run_Log(Base):
    """
    Run state transitions log.
    """

    __tablename__ = "run_logs"
    id = Column(Integer, primary_key=True)  # auto-increment
    run_id = Column(Integer, ForeignKey("runs.id"), nullable=False)
    status_from = Column(String(32), nullable=False)
    status_to = Column(String(32), nullable=False)
    timestamp = Column(DateTime, nullable=False, default=datetime.utcnow)
    reason = Column(String(1024), nullable=False, default="")
    sent = Column(Boolean, default=False, nullable=False)


class Task_Log(Base):
    __tablename__ = "task_logs"
    id = Column(Integer, primary_key=True)  # auto-increment
    cromwell_run_id = Column(String(36), nullable=False)
    execution_dir = Column(String(1024), nullable=False)
    status = Column(String(32), nullable=False)
    timestamp = Column(DateTime, nullable=False, default=datetime.utcnow)


class Transfer(Base):
    """
    Table of transfer tasks (sets of files to transfer).
    """

    __tablename__ = "transfers"
    id = Column(Integer, primary_key=True)
    status = Column(String(32), nullable=False)
    submitted = Column(DateTime, nullable=False, default=datetime.utcnow)
    updated = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    src_base_dir = Column(String(1024), nullable=False)
    dest_base_dir = Column(String(1024), nullable=False)
    manifest_json = Column(MEDIUMTEXT, nullable=False)
    reason = Column(String(256), nullable=True)


def create_all(engine, session):
    """Create all tables if not exist"""
    Base.metadata.create_all(engine)
    session.commit()
