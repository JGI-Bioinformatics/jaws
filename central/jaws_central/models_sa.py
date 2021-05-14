"""Vanilla Sqlalchemy ORM models, used by rpc_operations"""

from datetime.datetime import utcnow
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
    email = Column(String(64), nullable=False, unique=True)
    name = Column(String(64), nullable=True)
    is_admin = Column(Boolean, nullable=False, default=False)
    is_dashboard = Column(Boolean, nullable=False, default=False)
    jaws_token = Column(String(256), nullable=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<User {self.id}>"


class Run(Base):
    """Analysis runs are the execution of workflows on specific inputs."""

    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    user_id = Column(String(32), ForeignKey("users.id"), nullable=False)
    submission_id = Column(String(36), nullable=False)
    input_site_id = Column(String(8), nullable=False)
    compute_site_id = Column(String(8), nullable=False)
    wdl_file = Column(String(256), nullable=False)
    json_file = Column(String(256), nullable=False)
    tag = Column(String(256), nullable=True)
    status = Column(String(32), nullable=False)
    submitted = Column(DateTime, nullable=False, default=utcnow)
    updated = Column(
        DateTime, default=same_as("submitted"), onupdate=utcnow,
    )
    cromwell_run_id = Column(String(36), nullable=True)
    result = Column(String(32), nullable=True)
    upload_id = Column(Integer, nullable=True)
    download_id = Column(Integer, nullable=True)

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
    timestamp = Column(DateTime, nullable=False, default=utcnow)
    reason = Column(String(1024), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Run_Log {self.run_id}:{self.status_from}:{self.status_to}>"


class Xfer(Base):
    """Transfer tasks used by XferQueue class"""

    ___tablename__ = "xfers"
    id = Column(Integer, primary_key=True)
    user_id = Column(String(32), ForeignKey("users.id"), nullable=False)
    label = Column(String(32), nullable=True)
    src_endpoint = Column(String(36), nullable=False)
    dest_endpoint = Column(String(36), nullable=False)
    manifest = Column(Text, nullable=False)
    size_gb = Column(Float, nullable=False)
    priority_a = Column(Integer, nullable=False, default=1)
    priority_b = Column(Integer, nullable=False, default=1)
    status = column(String(16), nullable=False, default="created")
    submitted = column(DateTime, nullable=False, default=utcnow)
    updated = column(DateTime, nullable=False, default=utcnow)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Xfer {self.id}:{user_id}:{size_gb}GB:{status}>"
