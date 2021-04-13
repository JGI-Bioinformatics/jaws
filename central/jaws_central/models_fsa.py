"""Flask-SQLAlchemy db and models, used by Connexion/Flask servers."""

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
    email = db.Column(db.String(64), nullable=False)
    name = db.Column(db.String(64), nullable=True)
    is_admin = db.Column(db.Boolean, nullable=False, default=False)
    is_dashboard = db.Column(db.Boolean, nullable=False, default=False)
    jaws_token = db.Column(db.String(256), nullable=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

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
    user_id = db.Column(db.String(32), db.ForeignKey("users.id"), nullable=False)
    created = db.Column(db.DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = db.Column(
        db.DateTime,
        nullable=False,
        default=same_as("created"),
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
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<Workflow {self.name}:{self.version}>"


class Run(db.Model):
    """Analysis runs are the execution of workflows on specific inputs."""

    __tablename__ = "runs"
    id = db.Column(db.Integer, primary_key=True)
    submission_id = db.Column(db.String(36), nullable=False)
    cromwell_run_id = db.Column(db.String(36), nullable=True)
    result = db.Column(db.String(32), nullable=True)
    status = db.Column(db.String(32), nullable=False)
    user_id = db.Column(db.String(32), db.ForeignKey("users.id"), nullable=False)
    site_id = db.Column(db.String(8), nullable=False)
    submitted = db.Column(db.DateTime, nullable=False, default=datetime.datetime.utcnow)
    updated = db.Column(
        db.DateTime, default=same_as("submitted"), onupdate=datetime.datetime.utcnow,
    )
    input_site_id = db.Column(db.String(8), nullable=False)
    input_endpoint = db.Column(db.String(36), nullable=False)
    upload_task_id = db.Column(db.String(36), nullable=True)
    output_endpoint = db.Column(db.String(36), nullable=False)
    output_dir = db.Column(db.String(256), nullable=False)
    wdl_file = db.Column(db.String(256), nullable=False)
    json_file = db.Column(db.String(256), nullable=False)
    tag = db.Column(db.String(256), nullable=True)
    download_task_id = db.Column(db.String(36), nullable=True)

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
