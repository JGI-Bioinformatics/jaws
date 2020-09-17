"""
RPC operations for jaws-user service
"""

import logging
import sqlalchemy.exc
from datetime import datetime
import globus_sdk
from sqlalchemy.exc import SQLAlchemyError
from jaws_rpc.responses import success, failure
from jaws_user import config
from jaws_user.database import Session
from jaws_user.models import User

logger = logging.getLogger(__package__)


def get_user(params):
    """
    Retrieve a user's records by ID.

    :param params: contains "user_id"
    :type params: dict
    :return: valid JSON-RPC2 response
    :rtype: dict
    """
    user_id = params["user_id"]
    logger.debug(f"Get user {user_id}")
    session = Session()
    try:
        user = session.query(User).get(user_id)
    except SQLAlchemyError as e:
        return failure(500, f"Db error: {e}")
    if user is None:
        err = f"User {user_id} not found"
        logger.error(err)
        return failure(404, err)
    result = {
        "user_id": user.id,
        "email": user.email,
        "name": user.name,
        "is_admin": user.admin,
        "globus_id": user.globus_id,
        "auth_refresh_token": user.auth_refresh_token,
        "transfer_refresh_token": user.transfer_refresh_token,
    }
    session.close()
    return success(result)


def is_admin(params):
    """
    Check if a user has administrator access.

    :param params: contains "user_id"
    :type params: dict
    :return: JSON-RPC2 response with boolean result
    :rtype: dict
    """
    response = get_user(params)
    if "error" in response:
        return response
    else:
        return success(response["result"]["is_admin"])


def transfer_refresh_token(params):
    """
    Retrieve a user's Globus transfer refresh token.

    :param params: contains "user_id"
    :type params: dict
    :return: JSON-RPC2 response with str result
    :rtype: dict
    """
    response = get_user(params)
    if "error" in response:
        return response
    else:
        return success(response["result"]["transfer_refresh_token"])


def save_globus_tokens(user_id, auth_refresh_token, transfer_refresh_token):
    """Save a user's Globus tokens.  Only the refresh tokens are required with the RefreshTokenAuthorizer
    and the refresh tokens never expire.

    :param user_id: user ID
    :type user_id: str
    :param auth_refresh_token: Auth Refresh Token
    :type auth_refresh_token: str
    :param auth_refresh_token: Transfer Refresh Token
    :type transfer_refresh_token: str
    :return:
    """
    conf = config.conf

    # CHECK IF REGISTERED JAWS USER
    try:
        user_rec = db.session.query(User).get(user_id)
    except SQLAlchemyError as e:
        return (500, f"Db error: {e}")
    if user_rec is None:
        logger.error(f"No match for user {user_id}")
        return failure(401, "User db record not found")

    # UPDATE USER RECORD; ADD MISSING FIELDS IF NOT DEFINED
    if user_rec.globus_id is None:
        client = globus_sdk.NativeAppAuthClient(conf.get("GLOBUS", "client_id"))
        authorizer = globus_sdk.RefreshTokenAuthorizer(auth_refresh_token, client)
        auth_client = globus_sdk.AuthClient(authorizer=authorizer)
        user_info = auth_client.oauth2_userinfo()
        globus_id = user_info["sub"]
        name = user_info["name"]
        logger.debug(f"Defining name of user {user_id} as '{name}'")
        user_rec.name = name
        user_rec.globus_id = globus_id
        user_rec.auth_refresh_token = auth_refresh_token
        user_rec.transfer_refresh_token = transfer_refresh_token
    else:
        user_rec.auth_refresh_token = auth_refresh_token
        user_rec.transfer_refresh_token = transfer_refresh_token
    session.commit()
    return success()


## all RPC operations are defined in this dispatch table
operations = {
    "update_run_logs": {
        "function": update_run_logs,
        "required_parameters": [
            "site_id",
            "run_id",
            "status_from",
            "status_to",
            "timestamp",
        ],
    },
    "save_globus_tokens": {
        "function": save_globus_tokens,
        "required_parameters": [
            "user_id",
            "auth_refresh_token",
            "transfer_refresh_token",
        ],
    },
}
