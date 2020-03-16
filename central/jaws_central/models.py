"""
SQLAlchemy models for persistent data structures.
"""

import datetime
from sqlalchemy import Column, DateTime, String, Integer, Boolean, Text, ForeignKey, UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()


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
    """
    A registered user.
    """

    __tablename__ = "users"
    id = Column(String(36), primary_key=True)
    name = Column(String(64))
    email = Column(String(64), unique=True)
    is_admin = Column(
        Boolean, default=False, nullable=False
    )  # is True only if belong to jaws_admins Globus Group
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
    Once tagged as "released", it can't be edited or deleted, only deprecated.
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


class Run(Base):
    """
    Analysis runs are the execution of workflows on specific inputs.
    """

    __tablename__ = "runs"
    id = Column(Integer, primary_key=True)
    user_id = Column(String(36), ForeignKey("users.id"), nullable=False)
    site_id = Column(String(8), nullable=False)
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
#    workflow = relationship("Workflow", back_populates="runs")


def create_all(engine) -> None:
    """Create all tables which do not exist.
    """
    Base.metadata.create_all(engine)
