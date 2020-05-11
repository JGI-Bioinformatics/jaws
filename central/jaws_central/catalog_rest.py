"""
Workflows Catalog REST endpoints.
"""

from flask import abort, request
import markdown
import logging
from typing import Tuple
from jaws_central import catalog


logger = logging.getLogger(__package__)


def list_wdls(user: str) -> Tuple[dict, int]:
    """Retrieve workflows from database.

    :param user: Current user's ID
    :type user: str
    :return: Table of workflows
    :rtype: list
    """
    logger.debug("List workflows")
    try:
        result = catalog.list_wdls()
        return result, 200
    except catalog.CatalogDatabaseError as e:
        abort(500, f"Catalog error: {e}")


def get_versions(user: str, name: str) -> Tuple[dict, int]:
    """Returns a list of all (non-deprecated) versions of a particular workflow.

    :param user: Current user's ID
    :type user: str
    :param name: Name of the workflow
    :param str:
    :return: Table of all version of a workflow
    :rtype: list
    """
    logger.debug(f"Get version of workflow {name}")
    try:
        result = catalog.get_versions(name)
    except catalog.CatalogWorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}")
    except catalog.CatalogDatabaseError as e:
        abort(500, f"Catalog error: {e}")
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
    logger.debug(f"Get README of {name}:{version}")
    try:
        doc = catalog.get_doc(name, version)
    except catalog.CatalogWorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except catalog.CatalogDatabaseError as e:
        abort(500, f"Catalog error: {e}")
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
    logger.debug(f"Get WDL of {name}:{version}")
    try:
        wdl = catalog.get_wdl(name, version)
    except catalog.CatalogWorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except catalog.CatalogDatabaseError as e:
        abort(500, f"Catalog error: {e}")
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
    logger.debug(f"Release workflow, {name}:{version}")
    try:
        catalog.release_wdl(user, name, version)
    except catalog.CatalogWorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except catalog.CatalogAuthenticationError:
        abort(401, "Access denied; only the owner may release a workflow.")
    except catalog.CatalogDatabaseError as e:
        abort(500, f"Catalog error: {e}")
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
    logger.debug(f"Delete workflow, {name}:{version}")
    try:
        catalog.del_wdl(user, name, version)
    except catalog.CatalogWorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except catalog.CatalogAuthenticationError:
        abort(401, "Access denied; only the owner may delete a workflow.")
    except catalog.CatalogDatabaseError as e:
        abort(500, f"Catalog error: {e}")
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
    logger.debug(f"Update WDL of {name}:{version}")

    if "wdl_file" not in request.files or not request.files["wdl_file"]:
        abort(400, "Bad request: New WDL not provided")
    wdl_file = request.files["wdl_file"]
    new_wdl = wdl_file.read()
    if not new_wdl:
        abort(400, "Bad request: New WDL is empty")

    try:
        catalog.update_wdl(user, name, version, new_wdl)
    except catalog.CatalogWorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except catalog.CatalogAuthenticationError:
        abort(401, "Access denied; only the owner may update a workflow")
    except catalog.CatalogDatabaseError as e:
        abort(500, f"Catalog error: {e}")
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
    logger.debug(f"Update README of {name}:{version}")

    if "md_file" not in request.files or not request.files["md_file"]:
        abort(400, "Bad request: New MD not provided")
    md_file = request.files["md_file"]
    new_doc = md_file.read()
    if not new_doc:
        abort(400, "Bad request: New MD is empty")

    try:
        catalog.update_doc(user, name, version, new_doc)
    except catalog.CatalogWorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except catalog.CatalogAuthenticationError:
        abort(401, "Access denied; only the owner may update a workflow")
    except catalog.CatalogDatabaseError as e:
        abort(500, f"Catalog error: {e}")
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
    logger.debug(f"Add new workflow, {name}:{version}")

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
        catalog.add_wdl(user, name, version, new_wdl, new_doc)
    except catalog.CatalogDatabaseError as e:
        abort(500, f"Catalog error: {e}")
    return {"result": "OK", "name": name, "version": version}, 201
