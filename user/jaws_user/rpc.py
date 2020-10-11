"""
JAWS User Service RPC methods
"""

import logging
from jaws_rpc.responses import success, failure
from jaws_user.api import User, DatabaseError, UserNotFoundError

logger = logging.getLogger(__package__)


def get_user(params):
    """
    Retrieve a user's record.

    Required parameters: user_id
    Returns: dict of user information
    """
    try:
        user = User(params["user_id"])
    except Exception as error:
        return failure(error)
    return success(user)


def update_user(params):
    """
    Update user data.

    Required parameters: user_id, auth_refresh_token, transfer_refresh_token
    Return: None
    """
    try:
        user = User(params["user_id"])
        user.update_user(
            params["auth_refresh_token"],
            params["transfer_refresh_token"]
        )
    except Exception as error:
        return failure(error)
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
