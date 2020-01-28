"""
JAWS Workflows Catalog API
"""

import sys
import os
import json
import time
from datetime import datetime
from flask import make_response, abort, Flask, request, redirect, url_for
#from flask import current_app
import markdown
#from models import User, UserSchema, Workflow, WorkflowSchema
import models
from config import db

user_schema = models.UserSchema()
users_schema = models.UserSchema(many=True)
workflow_schema = models.WorkflowSchema()
workflows_schema = models.WorkflowSchema(many=True)

def list_wdls(user, release_filter=None):
    """
    Retrieve workflows from database.
    """
    workflows = models.Workflow.query.all() if release_filter is None else models.Workflow.query.filter_by(is_released=True).all()
    result = workflows_schema.dump(workflows)
    return { "workflows" : result }

def get_versions(user, name):
    """
    Returns a list of all (non-deprecated) versions of a particular workflow.
    """
    workflows = models.Workflow.query.filter_by(name=name, is_deprecated=False).all()
    result = workflows_schema.dump(workflows)
    return { "workflows" : result }

def get_doc(user, name, version):
    """
    Returns a workflow's README (stored in md format) as HTML.
    """
    workflow = models.Workflow.query.filter_by(name=name, version=version).one_or_none()
    if workflow is None:
        return {"message": "Workflow could not be found."}, 400
    return markdown.markdown(workflow.doc)

def get_wdl(user, name, version):
    """
    Returns a workflow's WDL
    """
    workflow = models.Workflow.query.filter_by(name=name, version=version).one_or_none()
    if workflow is None:
        return {"message": "Workflow could not be found."}, 400
    return workflow.wdl

def release_wdl(user, name, version):
    """
    Release a workflow (i.e. into production), which makes it's WDL immutable.
    """
    workflow = models.Workflow.query.filter_by(name=name, version=version).one_or_none()
    if workflow is None: abort(404, "Workflow not found")
    if workflow.uid != user: abort(401, "Access denied")
    if workflow.is_released is True: abort(400, "Workflow already released")
    workflow.is_released = True
    db.session.commit()
    return { "result" : "OK" }, 200

def del_wdl(user, name, version):
    """
    If a workflow has been released, it is marked as deprecated, otherwise it is purged from the db.
    """
    workflow = models.Workflow.query.filter_by(name=name, version=version).one_or_none()
    if workflow is None: abort(404, "Workflow not found")
    if workflow.uid != user: abort(401, "Access denied")
    if workflow.is_released is True:
        workflow.is_deprecated = True
    else:
        db.session.delete(workflow)
    db.session.commit()
    return { "result" : "OK" }, 200

def update_wdl(user, name, version):
    """
    Update the WDL file for a workflow.  This cannot be updated after release.
    """
    workflow = models.Workflow.query.filter_by(name=name, version=version).one_or_none()
    if workflow is None: abort(404, "Workflow not found")
    if workflow.uid != user: abort(401, "Access denied")
    if workflow.is_released is True: abort(400, "Released workflows cannot be updated")

    wdl_file = request.files['wdl_file']
    if not wdl_file: abort(400, "Bad request: WDL not provided")
    wdl = wdl_file.read()
    if not wdl: abort(400, "Bad request: WDL is empty")

    workflow.wdl = wdl
    db.session.commit()
    return { "result" : "OK" }, 200

def update_doc(user, name, version):
    """
    Update the doc file for a workflow.  This can be updated even after release.
    """
    workflow = models.Workflow.query.filter_by(name=name, version=version).one_or_none()
    if workflow is None: abort(404, "Workflow not found")
    if workflow.uid != user: abort(401, "Access denied")
    if workflow.is_released is True: abort(400, "Released workflows cannot be updated")

    md_file = request.files['md_file']
    if not md_file: abort(400, "Bad request: MD not provided")
    doc = md_file.read()
    if not doc: abort(400, "Bad request: MD is empty")

    workflow.doc = doc
    db.session.commit()
    return { "result" : "OK" }, 200

def add_wdl(user, name, version):
    """
    Add a new workflow to the catalog
    """
    # validate input
    name = name.lower().replace(" ", "_").replace(":", "__")
    version = version.lower().replace(" ", "_").replace(":", "__")
    workflow = models.Workflow.query.filter_by(name=name, version=version).one_or_none()
    if workflow is not None: abort(405, "Workflow already exists")
    wdl_file = request.files['wdl_file']
    if not wdl_file: abort(400, "Bad request: no WDL provided")
    wdl = wdl_file.read()
    if not wdl: abort(400, "Bad request: WDL is empty")
    md_file = request.files['md_file']
    if not md_file: abort(400, "Bad request: no README provided")
    doc = md_file.read()
    if not doc: abort(400, "Bad request: README is empty")
    now = datetime.now()

    # check if exists
    workflow = models.Workflow.query.filter_by(name=name, version=version).one_or_none()
    if workflow is not None: abort(400, "Workflow with that name:version already exists")

    # create new record
    workflow = models.Workflow(
        name=name,
        version=version,
        uid=user,
        date_created=now,
        is_released=False,
        is_deprecated=False,
        wdl=wdl,
        doc=doc 
    ) 
    db.session.add(workflow)
    db.session.commit()
    return { "result" : "OK" }, 201
