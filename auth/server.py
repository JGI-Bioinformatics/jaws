#!/usr/bin/env python

import sys
import os
import logging
import connexion
from connexion import request
from flask import abort
import requests
import configparser
import globus_sdk

import models

DEBUG = False
if "JAWS_DEBUG" in os.environ:
    DEBUG = os.environ["JAWS_DEBUG"]


# DB
db_session = None


# JAWS-AUTH CONFIG
if "JAWS_AUTH_CONFIG" not in os.environ:
    sys.exit("Env var $JAWS_AUTH_CONFIG not defined")
auth_config = configparser.ConfigParser()
auth_config.read_file(open(os.environ["JAWS_AUTH_CONFIG"]))

# GLOBUS
globus_client = globus_sdk.NativeAppAuthClient(auth_config["GLOBUS"]["client_id"])


def get_tokeninfo() -> dict:
    """
    OAuth2: validate token and return user info dictionary.  Abort wth code "401" on login failure.
    """
    # GET TOKEN
    try:
        # Discard "Bearer" keyword by splitting
        _, access_token = request.headers["Authorization"].split()
    except Exception:
        abort(401, "Authentication failure; invalid header")

    # UNPACK GLOBUS TOKENS
    if DEBUG:
        print("Got Globus super-token: %s" % (access_token,))
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
        db_session.query(models.User)
        .filter(models.User.auth_access_token == auth_access_token)
        .one_or_none()
    )
    if user is not None:
        # AUTH TOKEN MATCHES CACHED TOKEN, GET SCOPES; NO NEED TO QUERY GLOBUS
        scopes = ["jaws_users"]
        if user.is_admin:
            scopes.append("jaws_admins")
        if DEBUG:
            print("USER %s" % (user.id,))
        return {"uid": user.id, "scope": scopes}

    # CHECK IF VALID GLOBUS AUTH TOKEN
    if DEBUG:
        print("Checking if it's a valid Globus auth token")
    if not globus_client.oauth2_validate_token(auth_access_token)["active"]:
        abort(401, "Authentication failure")
    if DEBUG:
        print("Valid Globus auth token")

    # CHECK GLOBUS GROUPS TO SEE IF USER BELONGS TO "jaws_users"
    # THIS REQUIRES THE GLOBUS "GROUPS" TOKEN, NOT THE "AUTH" TOKEN
    scopes = []
    if DEBUG:
        print("Checking Globus groups")
    headers = {"Authorization": "Bearer %s" % (groups_access_token,)}
    r = requests.get(auth_config["GLOBUS"]["groups_url"], headers=headers)
    if r.status_code != requests.codes.ok:
        abort(r.status_code)
    result = r.json()
    if DEBUG:
        print(result)
    for group in result:
        group_name = group["name"]
        group_id = group["id"]
        if (
            group_name in auth_config["GROUPS"]
            and group_id == auth_config["GROUPS"][group_name]
        ):
            scopes.append(group_name)
    if "jaws_users" not in scopes:
        abort(401, "Not a member of jaws_user group")
    is_admin = True if "jaws_admins" in scopes else False

    # GET UID FROM GLOBUS INFO ENDPOINT USING "AUTH" TOKEN
    if DEBUG:
        print("Getting Globus user info")
    authorizer = globus_sdk.AccessTokenAuthorizer(auth_access_token)
    auth_client = globus_sdk.AuthClient(authorizer=authorizer)
    user_info = auth_client.oauth2_userinfo()
    id = user_info["sub"]
    name = user_info["name"]
    email = user_info["email"]
    if DEBUG:
        print("USER:\n\tid = %s\n\tname = %s\n\temail = %s" % (id, name, email))

    # CHECK IF USER EXISTS
    user = db_session.query(models.User).get(id)
    if user is None:
        # FIRST LOG-IN
        if DEBUG:
            print("New user; insert into db")
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
        db_session.add(new_user)
    else:
        # EXISTING USER; UPDATE TOKENS
        if DEBUG:
            print("Update existing user's tokens")
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
    db_session.commit()

    # DONE
    return {"uid": id, "scope": scopes}


def validate_globus_super_token(access_token):
    """
    Unpack set of Globus tokens and validate user vs Globus services.
    """


logging.basicConfig(level=logging.INFO)

app = connexion.FlaskApp(__name__)
app.add_api("swagger.yaml")

application = app.app
application.config.from_envvar("JAWS_DB_CONFIG")

db_session = models.init_db(application.config.get("SQLALCHEMY_DATABASE_URI"))


@application.teardown_appcontext
def shutdown_session(exception=None):
    db_session.remove()


if __name__ == "__main__":
    app.run(port=8080)
