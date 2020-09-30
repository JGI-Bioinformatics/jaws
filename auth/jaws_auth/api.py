import logging
from sqlalchemy.exc import SQLAlchemyError
from jaws_auth.db import db, User

logger = logging.getLogger(__package__)


def __select_user_by_token(access_token: str, session=None) -> dict:
    """
    Get user from db, if token exists.
    :param token: JAWS access token
    :type token: str
    :return: user record or None if user not found
    :rtype: dict
    """
    if not session:
        try:
            session = db.session()
        except SQLAlchemyError as error:
            logger.exception(f"Unable to get db session: {error}")
            raise

    # select row
    try:
        row = (
            session.query(User).filter(User.access_token == access_token).one_or_none()
        )
    except SQLAlchemyError as error:
        logger.exception(f"Error selecting user: {error}")
        raise
    return row


def get_user(access_token):
    """
    Get user record by token and return dict of uid and scopes.
    """
    try:
        row = __select_user_by_token(access_token)
    except SQLAlchemyError:
        raise
    if not row:
        return None
    user = {"uid": row.id, "scopes": ["user"]}
    if row.is_admin:
        user["scopes"].append("admin")
    return user
