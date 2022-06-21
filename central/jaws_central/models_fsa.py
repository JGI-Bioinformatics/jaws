"""
Flask-SQLAlchemy db and models, used by Connexion/Flask servers.
The tables are duplicated in the matching models.py file because
Flask-SqlAlchemy uses a different ORM base class than SqlAlchemy.
"""

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
    id = db.Column(db.String(32), primary_key=True)
    email = db.Column(db.String(64), nullable=True)
    name = db.Column(db.String(64), nullable=True)
    user_group = db.Column(db.String(16), nullable=True)
    is_admin = db.Column(db.Boolean, nullable=False, default=False)
    is_dashboard = db.Column(db.Boolean, nullable=False, default=False)
    jaws_token = db.Column(db.String(256), nullable=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<User {self.id}>"


class Transfer(db.Model):
    """
    Transfer tasks (i.e. uploads/downloads).
    If it's a Globus transfer, the globus_task_id will be specified.
    If it's an AWS-S3-copy transfer, the xfer_site_id will be specified.
    """

    __tablename__ = "transfers"
    id = db.Column(db.Integer, primary_key=True)
    status = db.Column(db.String(32), nullable=False, default="created")
    src_site_id = db.Column(db.String(8), nullable=False)
    src_base_dir = db.Column(db.String(128), nullable=False)
    dest_site_id = db.Column(db.String(8), nullable=False)
    dest_base_dir = db.Column(db.String(128), nullable=False)
    manifest_json = db.Column(db.String(64000), nullable=False)  # MEDIUMTEXT
    globus_transfer_id = db.Column(db.String(36), nullable=True)
    reason = db.Column(db.String(256), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Transfer {self.id}>"


class Run(db.Model):
    """Analysis runs are the execution of workflows on specific inputs."""

    __tablename__ = "runs"
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.String(32), db.ForeignKey("users.id"), nullable=False)
    submission_id = db.Column(db.String(36), nullable=False)
    max_ram_gb = db.Column(db.Integer, nullable=False, default=10)
    caching = db.Column(db.Boolean, nullable=False, default=True)
    input_site_id = db.Column(db.String(8), nullable=False)
    compute_site_id = db.Column(db.String(8), nullable=True)
    cromwell_run_id = db.Column(db.String(36), nullable=True)
    result = db.Column(db.String(32), nullable=True)
    status = db.Column(db.String(32), nullable=False, default="created")
    submitted = db.Column(db.DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = db.Column(
        db.DateTime,
        default=same_as("submitted"),
        onupdate=datetime.datetime.utcnow,
    )
    upload_id = db.Column(db.Integer, nullable=True)
    download_id = db.Column(db.Integer, nullable=True)
    wdl_file = db.Column(db.String(256), nullable=False)
    json_file = db.Column(db.String(256), nullable=False)
    tag = db.Column(db.String(256), nullable=True)
    manifest_json = db.Column(db.String(64000), nullable=False)  # MEDIUMTEXT
    webhook = db.Column(db.String(256), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Run {self.id}>"


class Run_Log(db.Model):
    """Run state transitions log"""

    __tablename__ = "run_logs"
    run_id = db.Column(db.Integer, db.ForeignKey("runs.id"), primary_key=True)
    status_from = db.Column(db.String(32), primary_key=True)
    status_to = db.Column(db.String(32), primary_key=True)
    timestamp = db.Column(db.DateTime, nullable=False, default=datetime.datetime.utcnow)
    reason = db.Column(db.String(1024), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Run_Log {self.run_id}:{self.status_from}:{self.status_to}>"
