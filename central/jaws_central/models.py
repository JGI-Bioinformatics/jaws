"""
Sqlalchemy ORM models
"""

import datetime
from sqlalchemy import (
    Column,
    DateTime,
    String,
    Integer,
    Boolean,
    ForeignKey,
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
    email = Column(String(64), nullable=True)
    name = Column(String(64), nullable=True)
    user_group = Column(String(16), nullable=True)
    is_admin = Column(Boolean, nullable=False, default=False)
    is_dashboard = Column(Boolean, nullable=False, default=False)
    jaws_token = Column(String(256), nullable=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<User {self.id}>"


class Transfer(Base):
    """
    Transfer tasks (i.e. uploads/downloads).
    If it's a Globus transfer, the globus_task_id will be specified.
    If it's an AWS-S3-copy transfer, the xfer_site_id will be specified.
    """

    __tablename__ = "transfers"
    id = Column(Integer, primary_key=True)
    status = Column(String(32), nullable=False, default="created")
    src_site_id = Column(String(8), nullable=False)
    src_base_dir = Column(String(128), nullable=False)
    dest_site_id = Column(String(8), nullable=False)
    dest_base_dir = Column(String(128), nullable=False)
    manifest_json = Column(String(64000), nullable=False)  # MEDIUMTEXT
    globus_transfer_id = Column(String(36), nullable=True)
    reason = Column(String(256), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Transfer {self.id}>"


class Run(Base):
    """Analysis runs are the execution of workflows on specific inputs."""

    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    user_id = Column(String(32), ForeignKey("users.id"), nullable=False)
    submission_id = Column(String(36), nullable=False)
    max_ram_gb = Column(Integer, nullable=False, default=10)
    caching = Column(Boolean, nullable=False, default=True)
    input_site_id = Column(String(8), nullable=False)
    compute_site_id = Column(String(8), nullable=True)
    cromwell_run_id = Column(String(36), nullable=True)
    result = Column(String(32), nullable=True)
    status = Column(String(32), nullable=False, default="created")
    submitted = Column(DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = Column(
        DateTime,
        default=same_as("submitted"),
        onupdate=datetime.datetime.utcnow,
    )
    upload_id = Column(Integer, nullable=True)
    download_id = Column(Integer, nullable=True)
    wdl_file = Column(String(256), nullable=False)
    json_file = Column(String(256), nullable=False)
    tag = Column(String(256), nullable=True)
    manifest_json = Column(String(64000), nullable=False)  # MEDIUMTEXT
    webhook = Column(String(256), nullable=True)

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


def create_all(engine, session):
    """Create all tables if not exist"""
    Base.metadata.create_all(engine)
    session.commit()
