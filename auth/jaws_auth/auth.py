#!/usr/bin/env python

import logging
import requests
import globus_sdk
from flask import abort, request
from jaws_auth import config, models

logger = logging.getLogger(__package__)
conf = config.JawsConfig()
db = conf.db
session = conf.session


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

    # UNPACK TOKENS
    (
        auth_access_token,
        auth_refresh_token,
        auth_expires_at_seconds,
        transfer_access_token,
        transfer_refresh_token,
        transfer_expires_at_seconds,
        groups_access_token,
        groups_refresh_token,
        groups_expires_at_seconds,
    ) = access_token.split(":")
    if not (auth_access_token and groups_access_token):
        abort(401, "Invalid super-token")

    # CHECK IF MATCHES TOKEN IN DB
    user = (
        session.query(models.User)
        .filter(models.User.auth_access_token == auth_access_token)
        .one_or_none()
    )
    if user is not None:
        # AUTH TOKEN MATCHES CACHED TOKEN, GET SCOPES; NO NEED TO QUERY GLOBUS
        scopes = ["jaws_users"]
        if user.is_admin:
            scopes.append("jaws_admins")
        return {"uid": user.id, "scope": scopes}

    # CHECK IF VALID GLOBUS AUTH TOKEN
    globus_client = globus_sdk.NativeAppAuthClient(conf.get_globus("client_id"))
    if not globus_client.oauth2_validate_token(auth_access_token)["active"]:
        abort(401, "Authentication failure")

    # CHECK GLOBUS GROUPS TO SEE IF USER BELONGS TO "jaws_users"
    # THIS REQUIRES THE GLOBUS "GROUPS" TOKEN, NOT THE "AUTH" TOKEN
    scopes = []
    headers = {"Authorization": "Bearer %s" % (groups_access_token,)}
    r = requests.get(conf.get_globus("groups_url"), headers=headers)
    if r.status_code != requests.codes.ok:
        abort(r.status_code)
    result = r.json()
    users_group_id = conf.get_globus("users_group")
    admins_group_id = conf.get_globus("admins_group")
    for group in result:
        if group["id"] == users_group_id or group["id"] == admins_group_id:
            scopes.append(group["name"])
    if "jaws_users" not in scopes:
        abort(401, "Not a member of jaws_user group")
    is_admin = True if "jaws_admins" in scopes else False

    # GET UID FROM GLOBUS INFO ENDPOINT USING "AUTH" TOKEN
    authorizer = globus_sdk.AccessTokenAuthorizer(auth_access_token)
    auth_client = globus_sdk.AuthClient(authorizer=authorizer)
    user_info = auth_client.oauth2_userinfo()
    id = user_info["sub"]
    name = user_info["name"]
    email = user_info["email"]

    # CHECK IF USER EXISTS
    user = session.query(models.User).get(id)
    if user is None:
        logger.info(f'Inserting new user, {name}, into db')
        new_user = models.User(
            id=id,
            name=name,
            email=email,
            is_admin=is_admin,
            auth_access_token=auth_access_token,
            auth_refresh_token=auth_refresh_token,
            auth_expires_at_seconds=auth_expires_at_seconds,
            transfer_access_token=transfer_access_token,
            transfer_refresh_token=transfer_refresh_token,
            transfer_expires_at_seconds=transfer_expires_at_seconds,
            groups_access_token=groups_access_token,
            groups_refresh_token=groups_refresh_token,
            groups_expires_at_seconds=groups_expires_at_seconds,
        )
        session.add(new_user)
    else:
        user.update(
            name,
            email,
            auth_access_token,
            auth_refresh_token,
            auth_expires_at_seconds,
            transfer_access_token,
            transfer_refresh_token,
            transfer_expires_at_seconds,
            groups_access_token,
            groups_refresh_token,
            groups_expires_at_seconds,
        )
    session.commit()

    return {"uid": id, "scope": scopes}
