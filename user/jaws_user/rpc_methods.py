"""
RPC operations for jaws-user service
"""

import logging
import globus_sdk
from sqlalchemy.exc import SQLAlchemyError
from jaws_rpc.responses import success, failure
from jaws_user import config
from jaws_user.db import Session, User

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


def update_user(user_id, auth_refresh_token, transfer_refresh_token):
    """Save a user's Globus tokens.  Only the refresh tokens are required with the RefreshTokenAuthorizer
    and the refresh tokens never expire.  This is used when new users grant the jaws app globus access.

    :param user_id: user ID
    :type user_id: str
    :param auth_refresh_token: Auth Refresh Token
    :type auth_refresh_token: str
    :param auth_refresh_token: Transfer Refresh Token
    :type transfer_refresh_token: str
    :return:
    """
    # GET USER INFO FROM GLOBUS
    client = globus_sdk.NativeAppAuthClient(config.conf.get("GLOBUS", "client_id"))
    authorizer = globus_sdk.RefreshTokenAuthorizer(auth_refresh_token, client)
    auth_client = globus_sdk.AuthClient(authorizer=authorizer)
    user_info = auth_client.oauth2_userinfo()
    globus_id = user_info["sub"]
    name = user_info["name"]
    email = user_info["email"]

    # INSERT OR UPDATE
    session = Session()
    try:
        user_rec = session.query(User).get(user_id)
    except SQLAlchemyError as e:
        return (500, f"Db error: {e}")
    if user_rec is None:
        # insert new
        new_user = User(
            name=name,
            email=email,
            globus_id=globus_id,
            auth_refresh_token=auth_refresh_token,
            transfer_refresh_token=transfer_refresh_token,
        )
        session.add(new_user)
    else:
        # update existing
        user_rec.name = name
        user_rec.globus_id = globus_id
        user_rec.auth_refresh_token = auth_refresh_token
        user_rec.transfer_refresh_token = transfer_refresh_token
    session.commit()
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
