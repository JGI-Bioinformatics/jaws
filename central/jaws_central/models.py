"""Database object and persistent object models."""

import datetime
from flask_sqlalchemy import SQLAlchemy


db = SQLAlchemy()


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


class User(db.Model):
    """Registered user"""

    __tablename__ = "users"
    id = db.Column(db.String(16), primary_key=True)
    email = db.Column(db.String(64), nullable=False)
    name = db.Column(db.String(64), nullable=True)
    is_admin = db.Column(db.Boolean, nullable=False, default=False)
    jaws_token = db.Column(db.String(256), nullable=False)
    globus_id = db.Column(db.String(36), nullable=True)
    auth_refresh_token = db.Column(db.String(256), nullable=True)
    transfer_refresh_token = db.Column(db.String(256), nullable=True)

    def __init__(self, *args, **kwargs):
        super(User, self).__init__(*args, **kwargs)

    def __repr__(self):
        return f"<User {self.id}>"


class Workflow(db.Model):
    """A workflow in the Catalog is comprised of WDL and MD files.
    Once tagged as "released", it can't be edited or deleted, only deprecated.
    """

    __tablename__ = "workflows"
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(32), nullable=False)
    version = db.Column(db.String(16), nullable=False, default="latest")
    user_id = db.Column(db.String(16), db.ForeignKey("users.id"), nullable=False)
    created = db.Column(db.DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = db.Column(
        db.DateTime,
        nullable=False,
        default=same_as("date_created"),
        onupdate=datetime.datetime.utcnow,
    )
    is_released = db.Column(db.Boolean, default=False, nullable=False)
    is_deprecated = db.Column(db.Boolean, default=False, nullable=False)
    wdl = db.Column(db.Text, nullable=False)
    doc = db.Column(db.Text, nullable=False)
    __table_args__ = (
        db.UniqueConstraint("name", "version", name="_workflow_name_version_uniq_cons"),
    )

    def __init__(self, *args, **kwargs):
        super(Workflow, self).__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Workflow {self.name}:{self.version}>"


class Run(db.Model):
    """Analysis runs are the execution of workflows on specific inputs."""

    __tablename__ = "runs"
    id = db.Column(db.Integer, primary_key=True)
    submission_id = db.Column(db.String(36), nullable=False)
    cromwell_id = db.Column(db.String(36), nullable=True)
    status = db.Column(db.String(16), nullable=False)
    user_id = db.Column(db.String(16), db.ForeignKey("users.id"), nullable=False)
    site_id = db.Column(db.String(8), nullable=False)
    submitted = db.Column(db.DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = db.Column(
        db.DateTime,
        default=same_as("date_submitted"),
        onupdate=datetime.datetime.utcnow,
    )
    input_site_id = db.Column(db.String(8), nullable=False)
    input_endpoint = db.Column(db.String(36), nullable=False)
    upload_task_id = db.Column(db.String(36), nullable=True)
    output_endpoint = db.Column(db.String(36), nullable=False)
    output_dir = db.Column(db.String(256), nullable=False)
    download_task_id = db.Column(db.String(36), nullable=True)

    def __init__(self, *args, **kwargs):
        super(Run, self).__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Run {self.id}>"
