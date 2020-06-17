"""Vanilla Sqlalchemy ORM models, used by rpc_operations"""

import datetime
from sqlalchemy import (
    Column,
    DateTime,
    String,
    Integer,
    Boolean,
    ForeignKey,
    Text,
    UniqueConstraint,
)
from jaws_central.database import Base


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


class User(Base):
    """Registered user"""

    __tablename__ = "users"
    id = Column(String(16), primary_key=True)
    email = Column(String(64), nullable=False)
    name = Column(String(64), nullable=True)
    is_admin = Column(Boolean, nullable=False, default=False)
    jaws_token = Column(String(256), nullable=False)
    globus_id = Column(String(36), nullable=True)
    auth_refresh_token = Column(String(256), nullable=True)
    transfer_refresh_token = Column(String(256), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<User {self.id}>"


class Workflow(Base):
    """A workflow in the Catalog is comprised of WDL and MD files.
    Once tagged as "released", it can't be edited or deleted, only deprecated.
    """

    __tablename__ = "workflows"
    id = Column(Integer, primary_key=True)
    name = Column(String(32), nullable=False)
    version = Column(String(16), nullable=False, default="latest")
    user_id = Column(String(16), ForeignKey("users.id"), nullable=False)
    created = Column(DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = Column(
        DateTime,
        nullable=False,
        default=same_as("date_created"),
        onupdate=datetime.datetime.utcnow,
    )
    is_released = Column(Boolean, default=False, nullable=False)
    is_deprecated = Column(Boolean, default=False, nullable=False)
    wdl = Column(Text, nullable=False)
    doc = Column(Text, nullable=False)
    __table_args__ = (
        UniqueConstraint("name", "version", name="_workflow_name_version_uniq_cons"),
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Workflow {self.name}:{self.version}>"


class Run(Base):
    """Analysis runs are the execution of workflows on specific inputs."""

    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    submission_id = Column(String(36), nullable=False)
    cromwell_id = Column(String(36), nullable=True)
    status = Column(String(16), nullable=False)
    user_id = Column(String(16), ForeignKey("users.id"), nullable=False)
    site_id = Column(String(8), nullable=False)
    submitted = Column(DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = Column(
        DateTime, default=same_as("date_submitted"), onupdate=datetime.datetime.utcnow,
    )
    input_site_id = Column(String(8), nullable=False)
    input_endpoint = Column(String(36), nullable=False)
    upload_task_id = Column(String(36), nullable=True)
    output_endpoint = Column(String(36), nullable=False)
    output_dir = Column(String(256), nullable=False)
    download_task_id = Column(String(36), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Run {self.id}>"


class Run_Log(Base):
    """Run state transitions log"""

    __tablename__ = "run_logs"
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey("runs.id"), nullable=False)
    status_from = Column(String(32), nullable=False)
    status_to = Column(String(32), nullable=False)
    timestamp = Column(DateTime, nullable=False, default=datetime.datetime.utcnow)
    reason = Column(String(1024), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Run_Log {self.id}>"


class Job_Log(Base):
    """A Run has many Tasks."""

    __tablename__ = "job_logs"
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey("runs.id"), nullable=False)
    cromwell_run_id = Column(String(36), nullable=False)
    cromwell_job_id = Column(Integer, nullable=False)
    task_name = Column(String(128), nullable=False)
    attempt = Column(Integer, nullable=False)
    status_from = Column(String(32), nullable=False)
    status_to = Column(String(32), nullable=False)
    timestamp = Column(DateTime, nullable=False)
    reason = Column(String(1024), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Job_Log {self.id}>"
