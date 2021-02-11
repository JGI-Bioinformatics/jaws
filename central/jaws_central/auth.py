import logging
import globus_sdk
from flask import abort, request
from sqlalchemy.exc import SQLAlchemyError
from jaws_central import config
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


def save_globus_tokens(
    user, auth_refresh_token, transfer_refresh_token,
):
    """Save a user's Globus tokens.  Only the refresh tokens are required with the RefreshTokenAuthorizer
    and the refresh tokens never expire.

    :param user: user ID
    :type user: str
    :param auth_refresh_token: Auth Refresh Token
    :type auth_refresh_token: str
    :param auth_refresh_token: Transfer Refresh Token
    :type transfer_refresh_token: str
    :return:
    """
    conf = config.conf

    # CHECK IF REGISTERED JAWS USER
    try:
        user_rec = db.session.query(User).get(user)
    except SQLAlchemyError as e:
        abort(500, f"Db error: {e}")
    if user_rec is None:
        logger.error(f"No match for user {user} in db")
        abort(401, "User db record not found")

    # UPDATE USER RECORD; ADD MISSING FIELDS IF NOT DEFINED
    if user_rec.globus_id is None:
        client = globus_sdk.NativeAppAuthClient(conf.get("GLOBUS", "client_id"))
        authorizer = globus_sdk.RefreshTokenAuthorizer(auth_refresh_token, client)
        auth_client = globus_sdk.AuthClient(authorizer=authorizer)
        user_info = auth_client.oauth2_userinfo()
        globus_id = user_info["sub"]
        name = user_info["name"]
        logger.debug(f"Defining name of user {user} as '{name}'")
        user_rec.name = name
        user_rec.globus_id = globus_id
        user_rec.auth_refresh_token = auth_refresh_token
        user_rec.transfer_refresh_token = transfer_refresh_token
    else:
        user_rec.auth_refresh_token = auth_refresh_token
        user_rec.transfer_refresh_token = transfer_refresh_token
    try:
        db.session.commit()
    except Exception as error:
        db.session.rollback()
        logger.exception(f"Error updating run with globus tokens: {error}")
        abort(500, f"Error saving Globus tokens; please try again later: {error}")


def _get_user_by_globus_id(globus_user_id):
    try:
        user = (
            db.session.query(User)
            .filter(User.globus_id == globus_user_id)
            .one_or_none()
        )
    except SQLAlchemyError as e:
        abort(500, f"Db error: {e}")
    if user is None:
        abort(404, "User not found")
    return user


def get_user_token(user, globus_user_id):
    """
    Return the JAWS token for the queried user' globus ID.
    This is restricted to 'dashboard' scope.

    :param user: user ID (i.e. dashboard user)
    :type user: str
    :param globus_user_id: The ID of the user logged into the dashboard
    :type globus_user_id: str
    :return: The user's JAWS access token.
    :rtype: str
    """
    query_user = _get_user_by_globus_id(globus_user_id)
    return {"jaws_token": query_user.jaws_token}
