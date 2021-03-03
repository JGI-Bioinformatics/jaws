import logging
from flask import abort, request
from sqlalchemy.exc import SQLAlchemyError
from jaws_central.models_fsa import db, User

logger = logging.getLogger(__package__)


def _get_bearer_token():
    """
    Extract the bearer token from the HTTP request header.
    """
    try:
        # Discard "Bearer" keyword by splitting
        _, access_token = request.headers["Authorization"].split()
    except Exception:
        abort(401, "Authentication failure; invalid header")
    return access_token


def _get_user_by_token(access_token):
    try:
        user = (
            db.session.query(User).filter(User.jaws_token == access_token).one_or_none()
        )
    except SQLAlchemyError as e:
        abort(500, f"Db error: {e}")
    if user is None:
        logger.info(f"Authentication failure; got token {access_token}")
        abort(401, "Authentication failure")
    return user


def get_tokeninfo() -> dict:
    """
    OAuth2: validate token and return user info dictionary.
    Abort on login failure.
    """
    access_token = _get_bearer_token()
    user = _get_user_by_token(access_token)

    logger.debug(f"User {user.id} token OK")
    scopes = ["user"]
    if user.is_admin:
        scopes.append("admin")
    if user.is_dashboard:
        scopes.append("dashboard")
    return {"uid": user.id, "scope": scopes}


def _get_user_by_email(email):
    try:
        user = (
            db.session.query(User)
            .filter(User.email == email)
            .one_or_none()
        )
    except SQLAlchemyError as e:
        abort(500, f"Db error: {e}")
    if user is None:
        abort(404, "User not found")
    return user


def get_user_token(user, email):
    """
    Return the JAWS token for the queried user' globus ID.
    This is restricted to 'dashboard' scope.

    :param user: user ID (i.e. dashboard user)
    :type user: str
    :param email: The email of the user logged into the dashboard
    :type email: str
    :return: The user's JAWS access token.
    :rtype: str
    """
    query_user = _get_user_by_email(email)
    return {"jaws_token": query_user.jaws_token}


def get_user(user):
    """
    Return current user's info.
    """
    try:
        user_rec = db.session.query(User).get(user)
    except SQLAlchemyError as e:
        abort(500, f"Db error: {e}")
    if user_rec is None:
        logger.error(f"No match for user {user} in db")
        abort(401, "User db record not found")
    result = {
        "uid": user,
        "name": user_rec.name,
        "email": user_rec.email,
    }
    return result
