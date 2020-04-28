"""
SQLAlchemy models for persistent data structures.
"""

import datetime
from sqlalchemy import (
    Column,
    DateTime,
    String,
    Integer,
    Boolean,
    Text,
    ForeignKey,
    UniqueConstraint,
)
from sqlalchemy.ext.declarative import declarative_base

# from sqlalchemy.orm import relationship

Base = declarative_base()


def same_as(column_name):
    """Function sets the default value of a column to the value in another column.

    :param column_name: Name of a column in a table.
    :type column_name: str
    :return: A function which returns the current value in another column.
    :rtype: function
    """

    def default_function(context):
        return context.current_parameters.get(column_name)

    return default_function


class User(Base):
    """Registered user"""

    __tablename__ = "users"
    id = Column(String(16), primary_key=True, unique=True)
    email = Column(String(64), nullable=False)
    name = Column(String(64), nullable=True)
    is_admin = Column(Boolean, nullable=False, default=False)
    jaws_token = Column(String(256), nullable=False)
    globus_id = Column(String(36), nullable=True)
    auth_refresh_token = Column(String(256), nullable=True)
    transfer_refresh_token = Column(String(256), nullable=True)


class Workflow(Base):
    """A workflow in the Catalog is comprised of WDL and MD files.
    Once marked as "released", a workflow cannot be changed or deleted, only deprecated.
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

    # ONE:MANY RELATIONSHIP
    # user = relationship("User", back_populates="workflows")


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
        DateTime, default=same_as("date_submitted"), onupdate=datetime.datetime.utcnow
    )
    wdl_file = Column(String(128), nullable=False)
    input_file = Column(String(128), nullable=False)
    input_site_id = Column(String(8), nullable=False)
    input_endpoint = Column(String(36), nullable=False)
    upload_task_id = Column(String(36), nullable=False)
    output_endpoint = Column(String(36), nullable=False)
    output_dir = Column(String(256), nullable=False)
    download_task_id = Column(String(36), nullable=True)

    # ONE:MANY RELATIONSHIPS


#    user = relationship("User", back_populates="runs")
#    workflow = relationship("Workflow", back_populates="runs")


def create_all(engine):
    """Create all tables if not exist"""
    Base.metadata.create_all(engine)
