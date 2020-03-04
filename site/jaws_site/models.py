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
from jaws_site import config

Base = declarative_base()


def same_as(column_name):
    """
    Function sets the default value of a column to the value in another column.
    """

    def default_function(context):
        return context.current_parameters.get(column_name)

    return default_function


class User(Base):
    """
    A registered user.
    """

    __tablename__ = "users"
    id = Column(String(36), primary_key=True)
    name = Column(String(64))
    email = Column(String(64), unique=True)
    is_admin = Column(
        Boolean, default=False, nullable=False
    )
    auth_access_token = Column(String(256))
    auth_refresh_token = Column(String(256))
    auth_expires_at_seconds = Column(Integer)
    transfer_access_token = Column(String(256))
    transfer_refresh_token = Column(String(256))
    transfer_expires_at_seconds = Column(Integer)
    groups_access_token = Column(String(256))
    groups_refresh_token = Column(String(256))
    groups_expires_at_seconds = Column(Integer)


class Workflow(Base):
    """
    A workflow in the Catalog is comprised of WDL and MD files.
    Once marked as "released", a workflow cannot be changed or deleted, only deprecated.
    """

    __tablename__ = "workflows"
    id = Column(Integer, primary_key=True)
    name = Column(String(32), nullable=False)
    version = Column(String(16), nullable=False, default="latest")
    user_id = Column(String(36), ForeignKey("users.id"), nullable=False)
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

    # MULTI-COLUMN CONSTRAINT
    __table_args__ = (
        UniqueConstraint("name", "version", name="_workflow_name_version_uniq_cons"),
    )

    # ONE:MANY RELATIONSHIP
    # user = relationship("User", back_populates="workflows")


class Site(Base):
    """
    Computing sites.
    """

    __tablename__ = "sites"
    id = Column(String(8), primary_key=True)
    endpoint = Column(String(36), nullable=False)
    basepath = Column(String(256), nullable=False)
    staging = Column(String(256), nullable=False)
    max_ram_gb = Column(Integer, nullable=False)
    max_transfer_gb = Column(Integer, nullable=True)
    has_dev_shm = Column(Boolean, nullable=False, default=False)


class Run(Base):
    """
    Analysis runs are the execution of workflows on specific inputs.
    """

    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    user_id = Column(String(36), ForeignKey("users.id"), nullable=False)
    site_id = Column(String(8), ForeignKey("sites.id"), nullable=False)
    submission_uuid = Column(String(36), nullable=False)
    status = Column(String(16), nullable=False)
    cromwell_id = Column(String(36), nullable=True)
    submitted = Column(DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = Column(
        DateTime, default=same_as("date_submitted"), onupdate=datetime.datetime.utcnow
    )
    upload_task_id = Column(String(36), nullable=False)
    download_task_id = Column(String(36), nullable=True)
    dest_endpoint = Column(String(36), nullable=False)
    dest_path = Column(String(256), nullable=False)

    # ONE:MANY RELATIONSHIPS
#    user = relationship("User", back_populates="runs")
#    site = relationship("Site", back_populates="runs")
#    workflow = relationship("Workflow", back_populates="runs")


def create_all():
    conf = config.JawsConfig()
    Base.metadata.create_all(conf.db)
