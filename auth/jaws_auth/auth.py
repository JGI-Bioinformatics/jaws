import logging
from flask import abort, request
from sqlalchemy.exc import SQLAlchemyError
from jaws_auth.db import db, User

logger = logging.getLogger(__package__)



def __select_user_by_token(token: str) -> dict:
    """
    Get user from db, if token exists.
    :param token: JAWS access token
    :type token: str
    :return: user record or None if user not found
    :rtype: dict
    """
    # get db session
    try:
        session = db.session()
    except Exception as error:
        logger.exception(f"Unable to get db session: {error}")
        raise
 
    # select row
    try:
        row = (
            db.session.query(User).filter(User.access_token == access_token).one_or_none()
        )
    except Exception as error:
        logger.exception(f"Error selecting user: {error}")
        raise
    return row

def __get_user_by_token(access_token):
    """
    Get user record by token and return dict of uid and scopes.
    """
    try:
        row = __select_user_by_token(access_token)
    except Exception as error:
        raise
    if not row:
        return None
    result = {
        "uid": row.id,
        "scopes": ["user"]
    }
    if row.is_admin:
        user["scopes"].append("admin")
    return result

def get_tokeninfo() -> dict:
    """
    OAuth2 REST endpoint: validate token in HTTP header and return user info dictionary.
    Abort on failure.
    """
    # Extract token from header.  Discard "Bearer" keyword by splitting
    try:
        _, access_token = request.headers["Authorization"].split()
    except Exception:
        abort(401, "Authentication failure; invalid header")

    # Get user record from db, if token exists
    try:
        user = __get_user_by_token(access_token)
    except Exception as error:
        abort(500, f"Auth db unavailable: {error}; please try again later")

    # Abort if invalid token or return uid and scopes
    if user is None:
        logger.info(f"Invalid token {access_token}")
        abort(401, "Authentication failure")
    else:
        return user
