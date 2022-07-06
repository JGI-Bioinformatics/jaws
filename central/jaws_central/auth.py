import logging
from flask import abort, request, current_app as app
from sqlalchemy.exc import SQLAlchemyError
import secrets
from jaws_central.models import User
import re
import requests
import json

logger = logging.getLogger(__package__)


OAUTH_TOKEN_LENGTH = 127


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
            app.session.query(User).filter(User.jaws_token == access_token).one_or_none()
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
        user = app.session.query(User).filter(User.email == email).one_or_none()
    except SQLAlchemyError as e:
        abort(500, f"Internal Database Error: {e}.")
    if user is None:
        abort(404, "User Not Found.")
    return user


def _get_json_from_sso(hash_code):
    base = "https://signon.jgi.doe.gov/api/app.sessions/"
    response = None
    try:
        response = requests.get(base + hash_code + ".json", allow_redirects=True)
    except Exception as err:
        logger.error(err)
        abort(response.status_code, f"{err}")
    try:
        sso_json = json.loads(response.content)
    except Exception as err:
        abort(403, f"SSO Hash was invalid: {err}")
    return sso_json


def _check_email_validity(email):
    regex = r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b"
    if re.fullmatch(regex, email):
        return True
    return False


def get_user_token(user, auth_key):
    """
    Return the JAWS token for the queried user' globus ID.
    This is restricted to 'dashboard' scope.

    :param user: user ID (i.e. dashboard user)
    :type user: str
    :param verify: either email id or sso token for authentication
    :type verify: str
    :return: The user's JAWS access token and/or sso json
    :rtype: dict
    """
    sso_json = {}
    if auth_key is None:
        abort(401, f"Authorization information missing: {auth_key}")
    if _check_email_validity(auth_key):
        email = auth_key
    elif not auth_key.isalnum():
        abort(401, f"Bad email address: {auth_key}.")
    else:
        sso_json = _get_json_from_sso(auth_key)
        if "user" not in sso_json or "email" not in sso_json["user"]:
            abort(401, f"No email address found in SSO JSON: {json.dumps(sso_json)}")
        email = sso_json["user"]["email"]
    query_user = _get_user_by_email(email)
    return {"jaws_token": query_user.jaws_token, "sso_json": sso_json}


def get_user(user):
    """
    Return current user's info.
    """
    try:
        user_rec = app.session.query(User).get(user)
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


def add_user(user) -> None:
    """
    Add user and return an OAuth2 token.
    """
    uid = request.form.get("uid", None).lower()
    name = request.form.get("name")
    user_group = request.form.get("user_group")
    email = request.form.get("email")
    admin = True if request.form.get("admin") == "True" else False

    a_user = app.session.query(User).get(uid)
    if a_user is not None:
        abort(400, "Cannot add user; uid already taken.")
    a_user = app.session.query(User).filter(User.email == email).one_or_none()
    if a_user is not None:
        abort(400, "Cannot add user; email already taken.")

    token = secrets.token_urlsafe(OAUTH_TOKEN_LENGTH)
    try:
        new_user = User(
            id=uid,
            name=name,
            user_group=user_group,
            jaws_token=token,
            email=email,
            is_admin=admin,
            is_dashboard=False,
        )
        app.session.add(new_user)
        app.session.commit()
        logger.info(f"{user} added new user {uid} ({email})")
    except Exception as error:
        app.session.rollback()
        logger.error(f"Failed to add user, {uid}: {error}")
        abort(500, f"Failed to add user, {uid}: {error}")
    return {"token": token}
