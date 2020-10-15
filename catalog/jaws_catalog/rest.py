"""
Workflows Catalog REST endpoints.


- The Workflows Catalog stores WDL files and accompanying markdown documentation.
- There may be multiple versions of a workflow of the same name.
- A workflow/version may be tagged as "released", in which case it cannot be deleted,
  but only tagged as "deprecated".
- Deprecated workflows do not appear in the catalog to users, but remain in the db.
- Workflows can only be updated/deleted by their owner.
"""

from flask import abort, request
import markdown
from typing import Tuple
from jaws_catalog.api import (
    Catalog,
    DatabaseError,
    PermissionsError,
    WorkflowNotFoundError,
    WorkflowImmutableError,
    InvalidInputError,
)
from jaws_catalog.database import db


def list_wdls(user: str) -> Tuple[dict, int]:
    """Retrieve workflows from database.

    :param user: Current user's ID
    :type user: str
    :return: Table of workflows
    :rtype: list
    """
    try:
        catalog = Catalog(db.session)
        result = catalog.list_wdls()
    except DatabaseError as error:
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
    try:
        catalog = Catalog(db.session)
        result = catalog.get_versions(name)
    except WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}")
    except DatabaseError as error:
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
    try:
        catalog = Catalog(db.session)
        doc = catalog.get_doc(name, version)
    except WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except DatabaseError as error:
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
    try:
        catalog = Catalog(db.session)
        wdl = catalog.get_wdl(name, version)
    except WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except DatabaseError as error:
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
    try:
        catalog = Catalog(db.session)
        catalog.release_wdl(user, name, version)
    except WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except PermissionsError:
        abort(401, "Access denied; only the owner may release a workflow.")
    except DatabaseError as error:
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
    try:
        catalog = Catalog(db.session)
        catalog.del_wdl(user, name, version)
    except WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except PermissionsError:
        abort(401, "Access denied; only the owner may delete a workflow.")
    except DatabaseError as error:
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
    # get WDL file
    if "wdl_file" not in request.files or not request.files["wdl_file"]:
        abort(400, "Bad request: New WDL not provided")
    wdl_file = request.files["wdl_file"]
    new_wdl = wdl_file.read()
    if not new_wdl:
        abort(400, "Bad request: New WDL is empty")

    # update entry
    try:
        catalog = Catalog(db.session)
        catalog.update_wdl(user, name, version, new_wdl)
    except WorkflowNotFoundError as error:
        abort(404, f"Workflow not found: {error}")
    except PermissionsError as error:
        abort(401, f"Access denied: {error}")
    except WorkflowImmutableError as error:
        abort(400, f"Action not allowed: {error}")
    except DatabaseError as error:
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
    # get markdown file
    if "md_file" not in request.files or not request.files["md_file"]:
        abort(400, "Bad request: New MD not provided")
    md_file = request.files["md_file"]
    new_doc = md_file.read()
    if not new_doc:
        abort(400, "Bad request: New MD is empty")

    # update entry
    try:
        catalog = Catalog(db.session)
        catalog.update_doc(user, name, version, new_doc)
    except WorkflowNotFoundError:
        abort(404, f"Workflow not found: {name}:{version}")
    except PermissionsError:
        abort(401, "Access denied; only the owner may update a workflow")
    except DatabaseError as error:
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
    # get WDL file
    if "wdl_file" not in request.files or not request.files["wdl_file"]:
        abort(400, "Bad request: New WDL not provided")
    wdl_file = request.files["wdl_file"]
    new_wdl = wdl_file.read()
    if not new_wdl:
        abort(400, "Bad request: New WDL is empty")

    # get markdown file
    if "md_file" not in request.files or not request.files["md_file"]:
        abort(400, "Bad request: New MD not provided")
    md_file = request.files["md_file"]
    new_doc = md_file.read()
    if not new_doc:
        abort(400, "Bad request: New MD is empty")

    # save new workflow
    try:
        catalog = Catalog(db.session)
        catalog.add_wdl(user, name, version, new_wdl, new_doc)
    except InvalidInputError as error:
        abort(400, f"Add workflow failed due to invalid input: {error}")
    except DatabaseError as error:
        abort(500, f"Catalog db error: {error}")
    return {"result": "OK"}, 201
