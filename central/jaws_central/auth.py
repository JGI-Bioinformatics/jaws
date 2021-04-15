import logging
from flask import abort, request
from sqlalchemy.exc import SQLAlchemyError
import secrets
from jaws_central.models_fsa import db, User
import requests
import json

logger = logging.getLogger(__package__)


OAUTH_TOKEN_LENGTH = 255


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
        user = db.session.query(User).filter(User.email == email).one_or_none()
    except SQLAlchemyError as e:
        abort(500, f"Db error: {e}")
    if user is None:
        abort(404, "User not found")
    return user


def _get_json_from_sso(hash_code):
    base = "https://signon.jgi.doe.gov/api/sessions/"
    response = None
    try:
        response = requests.get(base + hash_code + ".json", allow_redirects=True)
    except Exception as err:
        logger.error(err)
        abort(500, f"Error: {err}")
    if response is None:
        abort(401, "Response object is NoneType")
    return json.loads(response.content)


def get_user_token(user, hash_code):
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
    json_obj = _get_json_from_sso(hash_code)

    if "user" not in json_obj or "email" not in json_obj["user"]:
        abort(401, "No user or email information found in the JSON")

    query_user = _get_user_by_email(json_obj["user"]["email"])
    return {"jaws_token": query_user.jaws_token, "sso_json": json_obj}


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


def add_user(
    user: str, session, uid: str, name: str, email: str, admin: bool = False
) -> None:
    """
    Add user and return an OAuth2 token.
    """
    uid = uid.lower()
    a_user = db.session.query(User).get(uid)
    if a_user is not None:
        abort(400, f"Cannot add user {uid}; user.id already taken.")

    token = secrets.token_urlsafe(OAUTH_TOKEN_LENGTH)
    try:
        new_user = User(
            id=uid,
            name=name,
            jaws_token=token,
            email=email,
            is_admin=admin,
            is_dashboard=False,
        )
        db.session.add(new_user)
        db.session.commit()
        logger.info(f"Added new user {uid} ({email})")
    except Exception as error:
        db.session.rollback()
        logger.error(f"Failed to add user, {uid}: {error}")
        abort(500, f"Failed to add user, {uid}: {error}")
    return {"token": token}
