"""
Workflows Catalog REST endpoints.
"""

import os
from datetime import datetime
from flask import abort, request
import markdown
import logging
from sqlalchemy.exc import SQLAlchemyError
from jaws_central import config
from .models import Workflow

logger = logging.getLogger(__package__)
conf = config.JawsConfig()
session = conf.session


def list_wdls(user, release_filter=None):
    """Retrieve workflows from database.

    :param user: Current user's ID
    :type user: str
    :param release_filter: If flag set, then list only workflows tagged as "released"
    :type release_filter: bool, optional
    :return: Table of workflows
    :rtype: list
    """
    logger.info("List workflows")
    try:
        result = (
            session.query(Workflow)
            .filter(Workflow.is_released is True, Workflow.is_deprecated is False)
            .all()
            if release_filter
            else session.query(Workflow).filter(Workflow.is_deprecated is False).all()
        )
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Query failed")
    return result, 200


def get_versions(user, name):
    """Returns a list of all (non-deprecated) versions of a particular workflow.

    :param user: Current user's ID
    :type user: str
    :param name: Name of the workflow
    :param str:
    :return: Table of all version of a workflow
    :rtype: list
    """
    logger.info(f'Get version of workflow {name}')
    try:
        result = (
            session.query(Workflow)
            .filter(Workflow.name == name, Workflow.is_deprecated is False)
            .all()
        )
        return result, 200
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Query failed")


def get_doc(user, name, version):
    """Returns a workflow's README (stored in md format) as HTML.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return: The workflow's README document in markdown format
    :rtype: str
    """
    logger.info(f'Get README of {name}')
    try:
        workflow = (
            session.query(Workflow).filter_by(name=name, version=version).one_or_none()
        )
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Query failed")
    if workflow:
        return markdown.markdown(workflow.doc), 200
    else:
        abort(404, "Workflow not found")


def get_wdl(user, name, version):
    """Returns a workflow's WDL

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return: The workflow's specification document in WDL format
    :rtype: str
    """
    logger.info('Get WDL of {name}')
    try:
        workflow = (
            session.query(Workflow).filter_by(name=name, version=version).one_or_none()
        )
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Query failed")
    if workflow:
        return workflow.wdl, 200
    else:
        abort(404, "Workflow not found")


def release_wdl(user, name, version):
    """Tag a workflow as "released", which makes it's WDL immutable.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return: OK message upon success; abort otherwise
    :rtype: dict
    """
    logger.info(f'Release workflow, {name}:{version}')
    try:
        workflow = (
            session.query(Workflow).filter_by(name=name, version=version).one_or_none()
        )
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Query failed")
    if workflow is None:
        abort(404, "Workflow not found")
    if workflow.user != user:
        abort(401, "Access denied")
    if workflow.is_released:
        abort(400, "Workflow already released")
    workflow.is_released = True
    try:
        session.commit()
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Update failed")
    return {"result": "OK"}, 200


def del_wdl(user, name, version):
    """Delete a workflow.  If it was "released", then tags as "deprecated", rather than being purged from db.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return: OK message upon success; abort otherwise
    :rtype: dict
    """
    logger.info(f'Delete workflow, {name}:{version}')
    try:
        workflow = (
            session.query(Workflow).filter_by(name=name, version=version).one_or_none()
        )
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Query failed")
    if workflow is None:
        abort(404, "Workflow not found")
    if workflow.user != user:
        abort(401, "Access denied")
    if workflow.is_released:
        try:
            workflow.is_deprecated = True
            session.commit()
        except SQLAlchemyError as e:
            logger.error(e)
            abort(500, "Update failed")
    else:
        try:
            session.delete(workflow)
            session.commit()
        except SQLAlchemyError as e:
            logger.error(e)
            abort(500, "Delete failed")
    return {"result": "OK"}, 200


def update_wdl(user, name, version):
    """Update the WDL file for a workflow.  This cannot be updated after release.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :param wdl_file: new WDL document, from formData
    :type wdl_file: file
    :return: OK message upon success; abort otherwise
    :rtype: dict
    """
    logger.info(f'Update WDL of {name}:{version}')
    try:
        workflow = (
            session.query(Workflow).filter_by(name=name, version=version).one_or_none()
        )
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Query failed")
    if workflow is None:
        abort(404, "Workflow not found")
    if workflow.user != user:
        abort(401, "Access denied")
    if workflow.is_released:
        abort(400, "Workflow already released")

    wdl_file = request.files["wdl_file"]
    if not wdl_file:
        abort(400, "Bad request: WDL not provided")
    new_wdl = wdl_file.read()
    if not new_wdl:
        abort(400, "Bad request: WDL is empty")

    try:
        workflow.wdl = new_wdl
        session.commit()
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Update failed")
    os.remove(wdl_file)
    return {"result": "OK"}, 200


def update_doc(user, name, version):
    """Update the doc file for a workflow.  This can be updated even after release.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :param doc_file: new README document, from formData
    :type doc_file: file
    :return: OK message upon success; abort otherwise
    :rtype: dict
    """
    logger.info(f'Update README of {name}:{version}')
    try:
        workflow = (
            session.query(Workflow).filter_by(name=name, version=version).one_or_none()
        )
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Query failed")
    if workflow is None:
        abort(404, "Workflow not found")
    if workflow.user != user:
        abort(401, "Access denied")
    if workflow.is_released:
        abort(400, "Workflow already released")

    doc_file = request.files["doc_file"]
    if not doc_file:
        abort(400, "Bad request: README not provided")
    new_doc = doc_file.read()
    if not new_doc:
        abort(400, "Bad request: README is empty")

    try:
        workflow.doc = new_doc
        session.commit()
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Update failed")
    os.remove(doc_file)
    return {"result": "OK"}, 200


def add_wdl(user, name, version):
    """Add a new workflow to the catalog

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :param wdl_file: new WDL document, from formData
    :type wdl_file: file
    :param doc_file: new README document, from formData
    :type doc_file: file
    :return: OK message upon success; abort otherwise
    :rtype: dict
    """
    name = name.lower().replace(" ", "_").replace(":", "__")
    version = version.lower().replace(" ", "_").replace(":", "__")
    logger.info(f'Add new workflow, {name}:{version}')
    try:
        workflow = (
            session.query(Workflow).filter_by(name=name, version=version).one_or_none()
        )
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Query failed")
    if workflow is not None:
        abort(405, "Workflow already exists")
    wdl_file = request.files["wdl_file"]
    if not wdl_file:
        abort(400, "Bad request: no WDL provided")
    wdl = wdl_file.read()
    if not wdl:
        abort(400, "Bad request: WDL is empty")
    doc_file = request.files["doc_file"]
    if not doc_file:
        abort(400, "Bad request: no README provided")
    doc = doc_file.read()
    if not doc:
        abort(400, "Bad request: README is empty")
    now = datetime.now()
    workflow = Workflow(
        name=name,
        version=version,
        user=user,
        created=now,
        is_released=False,
        is_deprecated=False,
        wdl=wdl,
        doc=doc,
    )
    try:
        session.add(workflow)
        session.commit()
    except SQLAlchemyError as e:
        logger.error(e)
        abort(500, "Insert failed")
    return {"result": "OK"}, 201
