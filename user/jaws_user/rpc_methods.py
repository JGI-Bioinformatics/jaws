"""
JAWS User Service RPC methods
"""

import logging
from jaws_rpc.responses import success, failure
from jaws_user import api

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
    try:
        user = api.get_user(user_id)
    except api.DatabaseError as error:
        return failure(500, f"User service error: {error}")
    if user is None:
        err = f"User {user_id} not found"
        logger.error(err)
        return failure(404, err)
    return success(user)


def update_user(user_id, auth_refresh_token, transfer_refresh_token):
    """Save a user's Globus tokens.  Only the refresh tokens are required with the RefreshTokenAuthorizer
    and the refresh tokens never expire.  This is used when new users grant the jaws app globus access.
    If the user record exists, it will be updated, otherwise a new record shall be inserted.

    :param user_id: user ID
    :type user_id: str
    :param auth_refresh_token: Auth Refresh Token
    :type auth_refresh_token: str
    :param auth_refresh_token: Transfer Refresh Token
    :type transfer_refresh_token: str
    :return:
    """
    try:
        api.update_user(user_id, auth_refresh_token, transfer_refresh_token)
    except api.DatabaseError as error:
        return failure(500, f"User service error: {error}")
    return success()


# this dispatch table is used by jaws_rpc.rpc_server
rpc_methods = {
    "get_user": {"method": get_user, "required_parameters": ["user_id"]},
    "update_user": {
        "method": update_user,
        "required_parameters": [
            "user_id",
            "auth_refresh_token",
            "transfer_refresh_token",
        ],
    },
}
