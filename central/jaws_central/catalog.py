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
from models import Workflow
from config import db


def list_wdls(user, release_filter=None):
    """
    Retrieve workflows from database.
    """
    result = db.session.query(Workflow).filter(Workflow.is_released==True, Workflow.is_deprecated==False).all() if release_filter else db.session.query(Workflow).filter(Workflow.is_deprecated==False).all()
    return result


def get_versions(user, name):
    """
    Returns a list of all (non-deprecated) versions of a particular workflow.
    """
    result = db.session.query(Workflow).filter(name==name, is_deprecated==False).all()
    return result


def get_doc(user, name, version):
    """
    Returns a workflow's README (stored in md format) as HTML.
    """
    workflow = db.session.query(Workflow).filter_by(name=name, version=version).one_or_none()
    if workflow:
        return markdown.markdown(workflow.doc)
    else:
        abort(404, "Workflow not found")


def get_wdl(user, name, version):
    """
    Returns a workflow's WDL
    """
    workflow = db.session.query(Workflow).filter_by(name=name, version=version).one_or_none()
    if workflow:
        return workflow.wdl
    else:
        abort(404, "Workflow not found")
    workflow = models.Workflow.query.filter_by(name=name, version=version).one_or_none()


def release_wdl(user, name, version):
    """
    Release a workflow (i.e. into production), which makes it's WDL immutable.
    """
    workflow = db.session.query(Workflow).filter_by(name=name, version=version).one_or_none()
    if workflow is None: abort(404, "Workflow not found")
    if workflow.user != user: abort(401, "Access denied")
    if workflow.is_released: abort(400, "Workflow already released")
    workflow.is_released = True
    db.session.commit()
    return { "result" : "OK" }, 200


def del_wdl(user, name, version):
    """
    If a workflow has been released, it is marked as deprecated, otherwise it is purged from the db.
    """
    workflow = db.session.query(Workflow).filter_by(name=name, version=version).one_or_none()
    if workflow is None: abort(404, "Workflow not found")
    if workflow.user != user: abort(401, "Access denied")
    if workflow.is_released:
        workflow.is_deprecated = True
    else:
        db.session.delete(workflow)
    db.session.commit()
    return { "result" : "OK" }, 200


def update_wdl(user, name, version):
    """
    Update the WDL file for a workflow.  This cannot be updated after release.
    """
    workflow = db.session.query(Workflow).filter_by(name=name, version=version).one_or_none()
    if workflow is None: abort(404, "Workflow not found")
    if workflow.user != user: abort(401, "Access denied")
    if workflow.is_released: abort(400, "Workflow already released")

    wdl_file = request.files['wdl_file']
    if not wdl_file: abort(400, "Bad request: WDL not provided")
    new_wdl = wdl_file.read()
    if not new_wdl: abort(400, "Bad request: WDL is empty")

    workflow.wdl = new_wdl
    db.session.commit()
    os.remove(wdl_file)
    return { "result" : "OK" }, 200


def update_doc(user, name, version):
    """
    Update the doc file for a workflow.  This can be updated even after release.
    """
    workflow = db.session.query(Workflow).filter_by(name=name, version=version).one_or_none()
    if workflow is None: abort(404, "Workflow not found")
    if workflow.user != user: abort(401, "Access denied")
    if workflow.is_released: abort(400, "Workflow already released")

    doc_file = request.files['doc_file']
    if not doc_file: abort(400, "Bad request: README not provided")
    new_doc = doc_file.read()
    if not new_doc: abort(400, "Bad request: README is empty")

    workflow.doc = new_doc
    db.session.commit()
    os.remove(doc_file)
    return { "result" : "OK" }, 200


def add_wdl(user, name, version):
    """
    Add a new workflow to the catalog
    """
    # validate input
    name = name.lower().replace(" ", "_").replace(":", "__")
    version = version.lower().replace(" ", "_").replace(":", "__")
    workflow = db.session.query(Workflow).filter_by(name=name, version=version).one_or_none()
    if workflow is not None: abort(405, "Workflow already exists")
    wdl_file = request.files['wdl_file']
    if not wdl_file: abort(400, "Bad request: no WDL provided")
    wdl = wdl_file.read()
    if not wdl: abort(400, "Bad request: WDL is empty")
    doc_file = request.files['doc_file']
    if not doc_file: abort(400, "Bad request: no README provided")
    doc = doc_file.read()
    if not doc: abort(400, "Bad request: README is empty")
    now = datetime.now()

    # create new record
    workflow = Workflow(
        name=name,
        version=version,
        user=user,
        created=now,
        is_released=False,
        is_deprecated=False,
        wdl=wdl,
        doc=doc 
    ) 
    db.session.add(workflow)
    db.session.commit()
    return { "result" : "OK" }, 201
