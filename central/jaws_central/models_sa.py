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
    id = Column(String(32), primary_key=True)
    email = Column(String(64), nullable=False, unique=True)
    name = Column(String(64), nullable=True)
    is_admin = Column(Boolean, nullable=False, default=False)
    is_dashboard = Column(Boolean, nullable=False, default=False)
    jaws_token = Column(String(256), nullable=False)

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
    user_id = Column(String(32), ForeignKey("users.id"), nullable=False)
    created = Column(DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = Column(
        DateTime,
        nullable=False,
        default=same_as("created"),
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
    cromwell_run_id = Column(String(36), nullable=True)
    result = Column(String(32), nullable=True)
    status = Column(String(32), nullable=False)
    user_id = Column(String(32), ForeignKey("users.id"), nullable=False)
    site_id = Column(String(8), nullable=False)
    submitted = Column(DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = Column(
        DateTime, default=same_as("submitted"), onupdate=datetime.datetime.utcnow,
    )
    input_site_id = Column(String(8), nullable=False)
    input_endpoint = Column(String(36), nullable=False)
    upload_task_id = Column(String(36), nullable=True)
    output_endpoint = Column(String(36), nullable=False)
    output_dir = Column(String(256), nullable=False)
    wdl_file = Column(String(256), nullable=False)
    json_file = Column(String(256), nullable=False)
    tag = Column(String(256), nullable=True)
    download_task_id = Column(String(36), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Run {self.id}>"


class Run_Log(Base):
    """Run state transitions log"""

    __tablename__ = "run_logs"
    run_id = Column(Integer, ForeignKey("runs.id"), primary_key=True)
    status_from = Column(String(32), primary_key=True)
    status_to = Column(String(32), primary_key=True)
    timestamp = Column(DateTime, nullable=False, default=datetime.datetime.utcnow)
    reason = Column(String(1024), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Run_Log {self.run_id}:{self.status_from}:{self.status_to}>"
