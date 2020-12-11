"""
Workflows Catalog REST endpoints.
"""

from datetime import datetime
import logging
from sqlalchemy.exc import SQLAlchemyError
from typing import Tuple
from jaws_central.models_fsa import db, Workflow

logger = logging.getLogger(__package__)


class CatalogDatabaseError(Exception):
    def __init__(self, message):
        super().__init__(message)


class CatalogAuthenticationError(Exception):
    def __init__(self, message):
        super().__init__(message)


class CatalogWorkflowNotFoundError(Exception):
    def __init__(self, message):
        super().__init__(message)


class CatalogWorkflowImmutableError(Exception):
    def __init__(self, message):
        super().__init__(message)


class CatalogInvalidInputError(Exception):
    def __init__(self, message):
        super().__init__(message)


def list_wdls() -> Tuple[dict, int]:
    """Retrieve workflows from database.

    :return: Table of workflows
    :rtype: list
    """
    logger.debug("List workflows")
    try:
        query = db.session.query(Workflow).filter(Workflow.is_deprecated == 0).all()
    except SQLAlchemyError as error:
        logger.error(error)
        raise CatalogDatabaseError(error)
    result = []
    for row in query:
        is_released = "yes" if row.is_released else "no"
        result.append({
            "name": row.name,
            "version": row.version,
            "owner": row.user_id,
            "created": row.created,
            "last_updated": row.updated,
            "production_release": is_released})
    return result


def get_versions(name: str) -> Tuple[dict, int]:
    """Returns a list of all (non-deprecated) versions of a particular workflow.

    :param name: Name of the workflow
    :param str:
    :return: Table of all version of a workflow
    :rtype: list
    """
    logger.debug(f"Get version of workflow {name}")
    try:
        query = (
            db.session.query(Workflow)
            .filter(Workflow.name == name, Workflow.is_deprecated == 0)
            .all()
        )
    except SQLAlchemyError as error:
        logger.error(error)
        raise CatalogDatabaseError(error)
    result = {}
    for row in query:
        is_released = "yes" if row.is_released else "no"
        full_name = f"{row.name}:{row.version}"
        result[full_name] = {
            "name": row.name,
            "version": row.version,
            "owner": row.user_id,
            "created": row.created,
            "last_updated": row.updated,
            "production_release": is_released,
        }
    return result


def get_doc(name: str, version: str) -> Tuple[str, int]:
    """Returns a workflow's README (stored in md format) as HTML.

    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return: The workflow's README document in markdown format
    :rtype: str
    """
    logger.debug(f"Get README of {name}")
    try:
        workflow = (
            db.session.query(Workflow)
            .filter_by(name=name, version=version)
            .one_or_none()
        )
    except SQLAlchemyError as error:
        logger.error(error)
        raise CatalogDatabaseError(error)
    if workflow:
        return workflow.doc
    else:
        raise CatalogWorkflowNotFoundError("Workflow not found")


def get_wdl(name: str, version: str) -> Tuple[str, int]:
    """Returns a workflow's WDL

    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return: The workflow's specification document in WDL format
    :rtype: str
    """
    logger.debug("Get WDL of {name}")
    try:
        workflow = (
            db.session.query(Workflow)
            .filter_by(name=name, version=version)
            .one_or_none()
        )
    except SQLAlchemyError as error:
        logger.error(error)
        raise CatalogDatabaseError(error)
    if workflow:
        return workflow.wdl
    else:
        raise CatalogWorkflowNotFoundError("Workflow not found")


def release_wdl(user: str, name: str, version: str) -> Tuple[dict, int]:
    """Tag a workflow as "released", which makes it's WDL immutable.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return:
    """
    logger.debug(f"Release workflow, {name}:{version}")
    try:
        workflow = (
            db.session.query(Workflow)
            .filter_by(name=name, version=version)
            .one_or_none()
        )
    except SQLAlchemyError as error:
        logger.error(error)
        raise CatalogDatabaseError(error)
    if workflow is None:
        raise CatalogWorkflowNotFoundError("Workflow not found")
    if workflow.user_id != user:
        raise CatalogAuthenticationError("User is not owner")
    workflow.is_released = True
    try:
        db.session.commit()
    except SQLAlchemyError as error:
        db.session.rollback()
        logger.error(error)
        raise CatalogDatabaseError(error)


def del_wdl(user: str, name: str, version: str) -> Tuple[dict, int]:
    """Delete a workflow.  If it was "released", then tags as "deprecated", rather than being purged from db.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return:
    """
    logger.debug(f"Delete workflow, {name}:{version}")
    try:
        workflow = (
            db.session.query(Workflow)
            .filter_by(name=name, version=version)
            .one_or_none()
        )
    except SQLAlchemyError as error:
        logger.error(error)
        raise CatalogDatabaseError(error)
    if workflow is None:
        raise CatalogWorkflowNotFoundError("Workflow not found")
    if workflow.user_id != user:
        raise CatalogAuthenticationError("User is not owner")
    if workflow.is_released:
        try:
            workflow.is_deprecated = True
            db.session.commit()
        except SQLAlchemyError as error:
            db.session.rollback()
            logger.error(error)
            raise CatalogDatabaseError(error)
    else:
        try:
            db.session.delete(workflow)
            db.session.commit()
        except SQLAlchemyError as error:
            db.session.rollback()
            logger.error(error)
            raise CatalogDatabaseError(error)


def update_wdl(user: str, name: str, version: str, new_wdl: str) -> Tuple[dict, int]:
    """Update the WDL file for a workflow.  This cannot be updated after release.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :param new_wdl: new WDL document
    :type new_wdl: str
    :return: True if successful; False if not found
    :rtype: bool
    """
    logger.debug(f"Update WDL of {name}:{version}")
    if not name and version and new_wdl:
        raise CatalogInvalidInputError("Missing required field")
    name = name.lower().replace(" ", "_").replace(":", "__")
    version = version.lower().replace(" ", "_").replace(":", "__")
    try:
        workflow = (
            db.session.query(Workflow)
            .filter_by(name=name, version=version)
            .one_or_none()
        )
    except SQLAlchemyError as error:
        logger.error(error)
        raise CatalogDatabaseError(error)
    if workflow is None:
        raise CatalogWorkflowNotFoundError("Workflow not found, check the name/version")
    if workflow.user_id != user:
        raise CatalogAuthenticationError("You are not the owner of this workflow")
    if workflow.is_released is True:
        raise CatalogWorkflowImmutableError("The WDL of a 'released' workflow cannot be updated, only deleted")
    try:
        workflow.wdl = new_wdl
        db.session.commit()
    except SQLAlchemyError as error:
        db.session.rollback()
        logger.error(error)
        raise CatalogDatabaseError(error)


def update_doc(user: str, name: str, version: str, new_doc: str) -> Tuple[dict, int]:
    """Update the doc file for a workflow.  This can be updated even after release.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :param new_doc: New README in Markdown format
    :type new_doc: str
    :return:
    """
    logger.debug(f"Update README of {name}:{version}")
    if not name and version and new_doc:
        raise CatalogInvalidInputError("Missing required field")
    name = name.lower().replace(" ", "_").replace(":", "__")
    version = version.lower().replace(" ", "_").replace(":", "__")
    try:
        workflow = (
            db.session.query(Workflow)
            .filter_by(name=name, version=version)
            .one_or_none()
        )
    except SQLAlchemyError as error:
        logger.error(error)
        raise CatalogDatabaseError(error)
    if workflow is None:
        raise CatalogWorkflowNotFoundError("Workflow not found")
    if workflow.user_id != user:
        raise CatalogAuthenticationError("User is not owner")
    try:
        workflow.doc = new_doc
        db.session.commit()
    except SQLAlchemyError as error:
        db.session.rollback()
        logger.error(error)
        raise CatalogDatabaseError(error)


def add_wdl(
    user: str, name: str, version: str, new_wdl: str, new_doc: str
) -> Tuple[dict, int]:
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
    :return:
    """
    if not name and version and new_doc:
        raise CatalogInvalidInputError("Missing required field")
    name = name.lower().replace(" ", "_").replace(":", "__")
    version = version.lower().replace(" ", "_").replace(":", "__")
    logger.debug(f"Add new workflow, {name}:{version}")
    try:
        workflow = (
            db.session.query(Workflow)
            .filter_by(name=name, version=version)
            .one_or_none()
        )
    except SQLAlchemyError as error:
        logger.error(error)
        raise CatalogDatabaseError(error)
    if workflow is not None:
        raise CatalogInvalidInputError(
            "A workflow with that name:version already exists"
        )
    now = datetime.now()
    workflow = Workflow(
        name=name,
        version=version,
        user_id=user,
        created=now,
        updated=now,
        is_released=False,
        is_deprecated=False,
        wdl=new_wdl,
        doc=new_doc,
    )
    try:
        db.session.add(workflow)
        db.session.commit()
    except SQLAlchemyError as error:
        db.session.rollback()
        logger.error(error)
        raise CatalogDatabaseError(error)
