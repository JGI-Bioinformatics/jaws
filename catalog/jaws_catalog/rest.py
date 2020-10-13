"""
Workflows Catalog REST endpoints.
"""

from flask import abort, request
import markdown
import logging
from typing import Tuple
from jaws_catalog import api, db


logger = logging.getLogger(__package__)


def list_wdls(user: str) -> Tuple[dict, int]:
    """Retrieve workflows from database.

    :param user: Current user's ID
    :type user: str
    :return: Table of workflows
    :rtype: list
    """
    logger.info(f"User {user}: List workflows")
    try:
        catalog = api.Catalog(db.session)
        result = catalog.list_wdls()
    except api.DatabaseError as error:
        abort(500, f"Catalog db error: {error}")
    return result, 200


def get_versions(user: str, name: str) -> Tuple[dict, int]:
    """Returns a list of all (non-deprecated) versions of a particular workflow.

    :param user: Current user's ID
    :type user: str
    :param name: Name of the workflow
    :param str:
    :return: Table of all version of a workflow
    :rtype: list
    """
    logger.info(f"User {user}: Get version of workflow {name}")
    try:
        catalog = api.Catalog(db.session)
        result = catalog.get_versions(name)
    except api.WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}")
    except api.DatabaseError as error:
        abort(500, f"Catalog db error: {error}")
    return result, 200


def get_doc(user: str, name: str, version: str) -> Tuple[str, int]:
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
    logger.info(f"User {user}: Get README of {name}:{version}")
    try:
        catalog = api.Catalog(db.session)
        doc = catalog.get_doc(name, version)
    except api.WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except api.DatabaseError as error:
        abort(500, f"Catalog db error: {error}")
    return markdown.markdown(doc), 200


def get_wdl(user: str, name: str, version: str) -> Tuple[str, int]:
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
    logger.info(f"User {user}: Get WDL of {name}:{version}")
    try:
        catalog = api.Catalog(db.session)
        wdl = catalog.get_wdl(name, version)
    except api.WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except api.DatabaseError as error:
        abort(500, f"Catalog db error: {error}")
    return wdl, 200


def release_wdl(user: str, name: str, version: str) -> Tuple[dict, int]:
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
    logger.info(f"User {user}: Release workflow, {name}:{version}")
    try:
        catalog = api.Catalog(db.session)
        catalog.release_wdl(user, name, version)
    except api.WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except api.AuthenticationError:
        abort(401, "Access denied; only the owner may release a workflow.")
    except api.DatabaseError as error:
        abort(500, f"Catalog db error: {error}")
    return {"result": "OK"}, 200


def del_wdl(user: str, name: str, version: str) -> Tuple[dict, int]:
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
    logger.info(f"User {user}: Delete workflow, {name}:{version}")
    try:
        catalog = api.Catalog(db.session)
        catalog.del_wdl(user, name, version)
    except api.WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except api.AuthenticationError:
        logger.info(f"Not allowing user {user} to delete workflow {name}:{version}")
        abort(401, "Access denied; only the owner may delete a workflow.")
    except api.DatabaseError as error:
        abort(500, f"Catalog db error: {error}")
    return {"result": "OK"}, 200


def update_wdl(user: str, name: str, version: str) -> Tuple[dict, int]:
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
    logger.info(f"User {user}: Update WDL of {name}:{version}")

    if "wdl_file" not in request.files or not request.files["wdl_file"]:
        abort(400, "Bad request: New WDL not provided")
    wdl_file = request.files["wdl_file"]
    new_wdl = wdl_file.read()
    if not new_wdl:
        abort(400, "Bad request: New WDL is empty")

    try:
        catalog = api.Catalog(db.session)
        catalog.update_wdl(user, name, version, new_wdl)
    except api.WorkflowNotFoundError as error:
        abort(404, f"Workflow not found: {error}")
    except api.AuthenticationError as error:
        abort(401, f"Access denied: {error}")
    except api.WorkflowImmutableError as error:
        abort(400, f"Action not allowed: {error}")
    except api.DatabaseError as error:
        abort(500, f"Catalog db error: {error}")
    return {"result": "OK"}, 200


def update_doc(user: str, name: str, version: str) -> Tuple[dict, int]:
    """Update the doc file for a workflow.  This can be updated even after release.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :param md_file: new README document, from formData
    :type md_file: file
    :return: OK message upon success; abort otherwise
    :rtype: dict
    """
    logger.info(f"User {user}: Update README of {name}:{version}")

    if "md_file" not in request.files or not request.files["md_file"]:
        abort(400, "Bad request: New MD not provided")
    md_file = request.files["md_file"]
    new_doc = md_file.read()
    if not new_doc:
        abort(400, "Bad request: New MD is empty")

    try:
        catalog = api.Catalog(db.session)
        catalog.update_doc(user, name, version, new_doc)
    except api.WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except api.AuthenticationError:
        abort(401, "Access denied; only the owner may update a workflow")
    except api.DatabaseError as error:
        abort(500, f"Catalog db error: {error}")
    return {"result": "OK"}, 200


def add_wdl(user: str, name: str, version: str) -> Tuple[dict, int]:
    """Add a new workflow to the catalog

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :param wdl_file: new WDL document, from formData
    :type wdl_file: file
    :param md_file: new README document, from formData
    :type md_file: file
    :return: OK message upon success; abort otherwise
    :rtype: dict
    """
    logger.info(f"User {user}: Add new workflow, {name}:{version}")

    if "wdl_file" not in request.files or not request.files["wdl_file"]:
        abort(400, "Bad request: New WDL not provided")
    wdl_file = request.files["wdl_file"]
    new_wdl = wdl_file.read()
    if not new_wdl:
        abort(400, "Bad request: New WDL is empty")

    if "md_file" not in request.files or not request.files["md_file"]:
        abort(400, "Bad request: New MD not provided")
    md_file = request.files["md_file"]
    new_doc = md_file.read()
    if not new_doc:
        abort(400, "Bad request: New MD is empty")

    try:
        catalog = api.Catalog(db.session)
        catalog.add_wdl(user, name, version, new_wdl, new_doc)
    except api.InvalidInputError as error:
        abort(400, f"Add workflow failed due to invalid input: {error}")
    except api.DatabaseError as error:
        abort(500, f"Catalog db error: {error}")
    return {"result": "OK", "name": name, "version": version}, 201
