"""
JAWS User class provides Globus authentication and access to the supplementary JAWS RDB records, such as status of user's runs.
"""

import os
import grp
import stat
import sys
import getpass
import requests
import json
import time
import globus_sdk
from globus_sdk import AuthClient, AccessTokenAuthorizer

#from jaws_client import config

# INIT GLOBUS PARAMS
GLOBUS_CLIENT_ID = os.environ["GLOBUS_CLIENT_ID"]
globus_client = globus_sdk.NativeAppAuthClient(GLOBUS_CLIENT_ID)
globus = {
    'auth_url' : 'https://auth.globus.org/v2/oauth2/token',
    'user_info_url' : 'https://auth.globus.org/v2/oauth2/userinfo',
    'groups_url' : 'https://groups.api.globus.org/v2/groups/my_groups',
    'auth_service_name' : 'auth.globus.org',
    'transfer_service_name' : 'transfer.api.globus.org',
    'groups_service_name' : "04896e9e-b98e-437e-becd-8084b9e234a0" # NOTE NOT groups.api.globus.org !?!
}

class User:
    """
    User class
    """
    config_file = None
    config = None
    # NOTE: the following are not guaranteed to be defined (yet); use the corresponding getter methods instead
    _scopes = None
    _uid = None
    _info = None
    _name = None
    _username = None
    _auth_client = None
    _transfer_client = None
    AUTH = {} # tokens
    TRANSFER = {} # tokens
    GROUPS = {} # tokens

    def __init__(self, user_id):
        """
        Load user info from JAWS' RDB
        """
        assert(user_id)

        # LOAD USER INFO FROM JAWS DATABASE
        # select access_token where uuid = %s

        # AUTHENTICATE AND INIT USER OBJECT
        self.init()
        assert(self.auth_token)


    def _save_tokens(self, globus_service, globus_auth_data):
        """
        Save tokens in object, but do not write.
        """
        for token_name in ("access_token", "refresh_token", "expires_at_seconds"):
            self.globus_service[token_name] = str(globus_auth_data[token_name])


    def auth_client(self):
        """
        Get Globus auth client object.
        """
        if self._auth_client is not None: return self._auth_client
        access_token = self.config["AUTH"]["access_token"]
        # NOTE: not using automatic refresh; see: refresh_tokens()
        #authorizer = globus_sdk.RefreshTokenAuthorizer(refresh_token, globus_client, access_token=access_token, expires_at=expires_at_seconds)
        #authorizer = globus_sdk.RefreshTokenAuthorizer(refresh_token, globus_client)
        authorizer = globus_sdk.AccessTokenAuthorizer(access_token)
        self._auth_client = globus_sdk.AuthClient(authorizer=authorizer)
        return self._auth_client


    def transfer_client(self):
        """
        Get Globus transfer client object.
        """
        if self._transfer_client is not None: return self._transfer_client
        access_token = self.config["TRANSFER"]["access_token"]
        authorizer = globus_sdk.AccessTokenAuthorizer(access_token)
        self._transfer_client = globus_sdk.TransferClient(authorizer=authorizer)
        return self._transfer_client


    def auth_token(self):
        """
        Returns the Globus AUTH access token.  This will be passed along with each REST request to JAWS-Central and is sufficient to validate the user and lookup any needed information.
        """
        return self.config["AUTH"]["access_token"] 


    def scopes(self):
        """
        Get JAWS scopes from Globus Groups and verify user is authorized to use JAWS.
        """
        if self._scopes is not None: return self._scopes
        access_token = self.config["GROUPS"]["access_token"]
        if not access_token: sys.exit("Globus Groups access token not defined")
        headers={'Authorization': "Bearer %s" % (access_token,)}
        r = requests.get(globus['groups_url'], headers=headers)
        result = r.json()
        scopes = []
        for group in result:
            group_name = group["name"]
            # TODO use IDs not names
            #group_id = group["id"]
            #print("%s : %s" % (group_name, group_id)) # ECCE
            if group_name.startswith("jaws_"): scopes.append(group_name)
        self.scopes = scopes
        if "jaws_users" not in scopes: sys.exit("You are not authorized to use JAWS")


    def info(self):
        """
        Get user info from Globus.  The user's auxilary JAWS information is accessed separately (using SQLAlchemy).
        """
        if self._info is not None: return self._info
        self._info = self.auth_client().oauth2_userinfo()
        self._uid = self._info["sub"]
        self._name = self._info["name"]
        self._username = self._info["preferred_username"]
        return self._info


    def uid(self):
        """
        Returns the user's Globus "sub" (user UIID)
        """
        self.info()
        return self._uid


    def username(self):
        """
        Returns the user's Globus preferred_username
        """
        self.info()
        return self._username


    def name(self):
        """
        Returns the user's full name
        """
        self.info()
        return self._name


    def refresh_tokens(self):
        """
        Refresh any tokens nearing expiration.
        """
        epoch_time = int(time.time())
        any_changes = False
        for globus_service in ("AUTH", "TRANSFER", "GROUPS"):
            expiry = int(self.config[globus_service]["expires_at_seconds"])
            if expiry - epoch_time < 86400:
                any_changes = True
                refresh_token = self.config[globus_service]["refresh_token"]
                self._save_tokens(globus_service, globus_client.oauth2_refresh_token(refresh_token))
        if any_changes: self.config.write_new_config_file()

    def init(self):
        """
        Verify and refresh authorization.
        """
        self.refresh_tokens()
        self.scopes()
