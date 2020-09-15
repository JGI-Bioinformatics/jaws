"""Flask-SQLAlchemy db and models for the Workflows Catalog."""

import datetime
from flask_sqlalchemy import SQLAlchemy


db = SQLAlchemy()


class Workflow(db.Model):
    """A workflow in the Catalog is comprised of WDL and MD files.
    Once tagged as "released", it can't be edited or deleted, only deprecated.
    """

    __tablename__ = "workflows"
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(32), nullable=False)
    version = db.Column(db.String(16), nullable=False, default="latest")
    user_id = db.Column(db.String(32), nullable=False)
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
