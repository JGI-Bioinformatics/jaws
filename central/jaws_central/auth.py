import logging
import globus_sdk
from flask import abort, request
from sqlalchemy.exc import SQLAlchemyError
from jaws_central import config
from jaws_central.models_fsa import db, User

logger = logging.getLogger(__package__)


def get_tokeninfo() -> dict:
    """
    OAuth2: validate token and return user info dictionary.
    Abort on login failure.
    """
    try:
        # Discard "Bearer" keyword by splitting
        _, access_token = request.headers["Authorization"].split()
    except Exception:
        abort(401, "Authentication failure; invalid header")

    try:
        user = (
            db.session.query(User).filter(User.jaws_token == access_token).one_or_none()
        )
    except SQLAlchemyError as e:
        abort(500, f"Db error: {e}")
    if user is None:
        logger.info(f"Authentication failure; got token {access_token}")
        abort(401, "Authentication failure")
    else:
        logger.debug(f"User {user.id} token OK")
        scopes = ["user"]
        if user.is_admin:
            scopes.append("admin")
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
    db.session.commit()
