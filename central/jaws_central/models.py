#!/usr/bin/env python

"""
JAWS persistent data structures: SQLAlchemy ORM Models and Mashmallow Schemas
"""

import datetime
import pytz
import enum
from marshmallow import fields
from marshmallow_enum import EnumField
from config import db, ma

#from sqlalchemy import Column, DateTime, String, Integer, Boolean, Text, create_engine, ForeignKey
#from sqlalchemy.ext.declarative import declarative_base
#from sqlalchemy.orm import scoped_session, sessionmaker
#Base = declarative_base()


####################
## SQLALCHEMY MODELS

def same_as(column_name):
    """
    Function sets the default value of a column to the value in another column
    """
    def default_function(context):
        return context.current_parameters.get(column_name)
    return default_function


class User(db.Model):
    """
    A registered user may have many Workflows and Runs.
    """
    __tablename__ = 'users'
    uid = db.Column(db.Integer, primary_key=True) # globus uuid
    name = db.Column(db.String(32), nullable=False)
    email = db.Column(db.String(64), nullable=False, unique=True)
    is_active = db.Column(db.Boolean, default=True)
    date_created = db.Column(db.DateTime, nullable=False, default=datetime.datetime.utcnow)
    date_expires = db.Column(db.DateTime, nullable=True, default=None)
    time_zone = db.Column(db.String(36), nullable=False)
    auth_token = db.Column(db.String(256), nullable=False)
    auth_scope = db.Column(db.String(20), nullable=True)
    globus_endpoint = db.Column(db.String(36), nullable=False)
    globus_path = db.Column(db.String(256), nullable=False)


class Workflow(db.Model):
    """
    A workflow is comprised of WDL and MD files.  A workflow cannot be changed or deleted after it has been released, only deprecated.
    """
    __tablename__ = 'workflows'
    workflow_id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(32), nullable=False)
    version = db.Column(db.String(16), nullable=False, default="latest")
    uid = db.Column(db.Integer, db.ForeignKey("users.uid"), nullable=False)
    date_created = db.Column(db.DateTime, nullable=False, default=datetime.datetime.utcnow)
    date_updated = db.Column(db.DateTime, nullable=False, default=same_as("date_created"), onupdate=datetime.datetime.utcnow)
    is_released = db.Column(db.Boolean, default=False)
    is_deprecated = db.Column(db.Boolean, default=False)
    wdl = db.Column(db.Text, nullable=False)
    doc = db.Column(db.Text, nullable=False)

    #user = db.relationship("User", backref="workflows")
    #user = db.relationship("User", backref="workflows", lazy="dynamic")
    #user = db.relationship("User", back_populates="workflows")

    __table_args__ = (db.UniqueConstraint('name', 'version', name='_workflow_name_version_uniq_cons'),)

#User.workflows = db.relationship("Workflow", order_by=Workflow.workflow_id, back_populates = "user")
#User.workflows = db.relationship("Workflow", backref="workflows")



class Site(db.Model):
    """
    Computing sites, with some info on available resources.
    """
    __tablename__ = "sites"
    site_id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(8), nullable=False)
    globus_endpoint = db.Column(db.String(36), nullable=False)
    globus_path = db.Column(db.String(256), nullable=False)
    max_ram_mb = db.Column(db.Integer)
    has_dev_shm = db.Column(db.Boolean)


class States(enum.Enum):
    """
    This enum describes the allowed scopes.
    """
    submitted = 1
    uploading = 2
    preprocessing = 3
    queued = 4
    running = 5
    postprocessing = 6
    downloading = 7
    completed = 8
    aborted = 99
    canceled = 100


class Run(db.Model):
    """
    Runs
    """
    __tablename__ = "runs"
    run_id = db.Column(db.Integer, primary_key=True)
    uid = db.Column(db.Integer, db.ForeignKey("users.uid"))
    site_id = db.Column(db.Integer, db.ForeignKey("sites.site_id"))
    analysis_id = db.Column(db.String(36), nullable=False)
    date_submitted = db.Column(db.DateTime(), nullable=False, default=datetime.datetime.utcnow)
    date_updated = db.Column(db.DateTime, default=same_as("date_submitted"), onupdate=datetime.datetime.utcnow)
    state = db.Column(db.Enum(States), nullable=False)
    #user = db.relationship("User", backref="runs", lazy="dynamic")



####################################
## MARSHMALLOW SERIALIZATION SCHEMAS

def must_not_be_blank(data):
    """
    Custom validator
    """
    if not data:
        raise ValidationError("Data not provided.")

class UserSchema(ma.ModelSchema):
    """
    Marshmallow schema for User class.
    """
    uid = fields.Int(dump_only=True)
    name = fields.Str(dump_only=True)
    email = fields.Str(dump_only=True)
    is_active = fields.Boolean(dump_only=True)
    date_created = fields.DateTime(dump_only=True)
    date_expires = fields.DateTime(dump_only=True)
    time_zone = fields.Str(dump_only=True)
    auth_token = fields.Str(dump_only=True)
    auth_scope = fields.Str(dump_only=True)
    globus_endpoint = fields.Str(dump_only=True)
    globus_path = fields.Str(dump_only=True)

    #workflows = fields.Nested("UserWorkflowSchema", default=[], many=True)
    #TODO runs = fields.Nested("UserRunSchema", default=[], many=True)

    class Meta:
        model = User
        sqla_session = db.session
        dateformat = '%Y-%m-%dT%H:%M:%S+0:00'


class WorkflowSchema(ma.ModelSchema):
    """
    Marshmallow schema for Workflow
    """
    workflow_id = fields.Int(dump_only=True)
    name = fields.Str(required=True, validate=must_not_be_blank)
    version = fields.Str(required=True, validate=must_not_be_blank)
    uid = fields.Int()
    date_created = fields.DateTime(dump_only=True)
    date_updated = fields.DateTime()
    is_released = fields.Boolean()
    is_deprecated = fields.Boolean()
    wdl = fields.String()
    doc = fields.String()
    user = fields.Nested("WorkflowPersonSchema", default=None)

    class Meta:
        model = Workflow
        sqla_session = db.session

class SiteSchema(ma.ModelSchema):
    """
    Marshmallow Site
    """
    site_id = fields.Int(dump_only=True)
    name = fields.Str(dump_only=True)
    globus_endpoint = fields.Str(dump_only=True)
    globus_path = fields.Str(dump_only=True)
    max_ram_mb = fields.Int(dump_only=True)
    has_dev_shm = fields.Boolean(dump_only=True)

    class Meta:
        model = Site
        sqla_session = db.session

    #runs = fields.Nested("SiteRunSchema", default=[], many=True)

class RunSchema(ma.ModelSchema):
    run_id = fields.Int()
    uid = fields.Int()
    site_id = fields.Int()
    analysis_id = fields.Str()
    date_submitted = fields.DateTime()
    date_updated = fields.DateTime()
    state = EnumField(States)

    class Meta:
        model = Run
        sqla_session = db.session

    #user = fields.Nested("RunPersonSchema", default=None)

























#class UserRunSchema(ma.ModelSchema):
#    """
#    Marshmallow schema of a user's runs.
#    """
#    def __init__(self, **kwargs):
#        super().__init__(strict=True, **kwargs)
#
#    # run fields:
#    run_id = fields.Int()
#    site_name = fields.nested(SiteSchema, only(name))
#    analysis_id = fields.Str()
#    submitted = fields.Str()
#    updated = fields.Str()
#    state = EnumField(States)

#class UserWorkflowSchema(ma.ModelSchema):
#    """
#    Marshmallow schema of a user's workflows.
#    """
#    def __init__(self, **kwargs):
#        super().__init__(strict=True, **kwargs)
#
#    # workflow fields:
#    workflow_id = fields.Int()
#    name = fields.Str()
#    version = fields.Str()
#    created = fields.Str()
#    updated = fields.Str()
#    released = fields.Boolean()
#    deprecated = fields.Boolean()

#class RunUserSchema(ma.ModelSchema):
#    """
#    Marshmallow schema of a run's user.
#    """
#    def __init__(self, **kwargs):
#        super().__init__(strict=True, **kwargs)
#
#    # user fields:
#    uid = fields.Int()
#    uname = fields.Str()
#    fname = fields.Str()
#    lname = fields.Str()
#    email = fields.Str()
#    time_zone = fields.Str()

#class WorkflowUserSchema(ma.ModelSchema):
#    """
#    Marshmallow schema of a workflow's user/owner.
#    """
#    def __init__(self, **kwargs):
#        super().__init__(strict=True, **kwargs)
#
#    # user fields:
#    uid = fields.Int()
#    uname = fields.Str()
#    fname = fields.Str()
#    lname = fields.Str()
#    email = fields.Str()
#    time_zone = fields.Str()

#class SiteRunSchema(ma.ModelSchema):
#    """
#    Marshmallow schema for a site's runs.
#    """
#    def __init__(self, **kwargs):
#        super().__init__(strict=True, **kwargs)
#
#    # run fields:
#    run_id = fields.Int()
#    uid = fields.Int()
#    analysis_id = fields.Str()
#    submitted = fields.Str()
#    updated = fields.Str()
#    state = EnumField(States)

####

#def init_db(uri):
#    engine = create_engine(uri, convert_unicode=True)
#    db_session = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=engine))
#    Base.query = db_session.query_property()
#    Base.metadata.create_all(bind=engine)
#    return db_session

if __name__ == "__main__":
    """
    Build database
    """
    db.create_all()
