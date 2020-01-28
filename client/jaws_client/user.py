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

from jaws_client import config

DEBUG = False
if "JAWS_DEBUG" in os.environ: DEBUG = True

# INIT GLOBUS PARAMS
JAWS_BRANCH = os.environ["JAWS_BRANCH"]
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
    _email = None
    _auth_client = None
    _transfer_client = None

    def __init__(self):
        """
        Load params from config file
        """
        # LOAD CONFIG
        if "HOME" not in os.environ: sys.exit('Env var "HOME" not defined')
        self.config_file = os.path.join(os.environ["HOME"], ".jaws.%s.ini" % (JAWS_BRANCH,))
        self.config = config.Config(path=self.config_file, template="config.ini")

        # AUTHENTICATE AND INIT USER OBJECT
        auth_at = self.config["AUTH"]["access_token"]
        if not auth_at:
            print("Welcome, new user")
            if input("Do you have a Globus account? [y/n]: ").upper()[0] != "Y": sys.exit("You must create a Globus login first.")
            if input("Have you been added to the JAWS group? [y/n]: ").upper()[0] != "Y": sys.exit("You must email an admin with your Globus username and request access to JAWS.")
            self.login()
        elif not globus_client.oauth2_validate_token(auth_at)['active']:
            print("Your authentication has expired.")
            self.login()
        else:
            self.init()
        assert(self.auth_token)


    def login(self):
        """
        Authenticate using Globus OAuth2 and init object params.  Sets self.auth_token upon success; exits otherwise.
        """
        # GLOBUS LOGIN
        requested_scopes = "openid profile email urn:globus:auth:scope:transfer.api.globus.org:all urn:globus:auth:scope:groups.api.globus.org:view_my_groups_and_memberships"
        globus_client.oauth2_start_flow(requested_scopes=requested_scopes, refresh_tokens=True)
        authorize_url = globus_client.oauth2_get_authorize_url()
        print("Please go to this URL and log in\n(HINT: try command-click on the link):\n\n%s\n" % (authorize_url,))
        auth_code = input('Then paste the authorization code here: ').strip()
        try:
            token_response = globus_client.oauth2_exchange_code_for_tokens(auth_code)
        except:
            sys.exit("Authentication failed")

        # SAVE TOKENS
        self._save_tokens_login("AUTH", token_response.by_resource_server[globus['auth_service_name']])
        self._save_tokens_login("TRANSFER", token_response.by_resource_server[globus['transfer_service_name']])
        self._save_tokens_login("GROUPS", token_response.by_resource_server[globus['groups_service_name']])

        # SET JAWS SCOPES FROM GLOBUS GROUPS; EXITS IF NOT AN AUTHORIZED JAWS USER.
        self.scopes()

        # WRITE NEW TOKENS TO USER CONFIG FILE
        self.config.write_new_config_file()
        print("Welcome, %s" % (self.name(),))


    def _save_tokens_login(self, globus_service, globus_auth_data):
        """
        Save tokens in config object, but do not write.
        """
        if DEBUG: print("Saving Globus %s token" % (globus_service,))
        self.config[globus_service]["access_token"] = globus_auth_data["access_token"]
        self.config[globus_service]["refresh_token"] = globus_auth_data["refresh_token"]
        self.config[globus_service]["expires_at_seconds"] = str(globus_auth_data["expires_at_seconds"])


    def _save_tokens_refresh(self, globus_service, globus_auth_data):
        """
        Save tokens in config object, but do not write.
        """
        if DEBUG: print("Saving Globus %s token" % (globus_service,))
        self.config[globus_service]["access_token"] = globus_auth_data["access_token"]
        self.config[globus_service]["refresh_token"] = globus_auth_data["refresh_token"]
        self.config[globus_service]["expires_at_seconds"] = str(int(time.time() + globus_auth_data["expires_in"]))


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
        print(self._info) # ECCE
        self._uid = self._info["sub"]
        self._name = self._info["name"]
        self._email = self._info["email"]
        return self._info


    def uid(self):
        """
        Returns the user's Globus "sub" (user UIID)
        """
        self.info()
        return self._uid


    def email(self):
        """
        Returns the user's email
        """
        self.info()
        return self._email


    def name(self):
        """
        Returns the user's full name
        """
        self.info()
        return self._name

    def header(self):
        """
        Return a dict to be used as HTTP OAuth2 header containing the current user's authentication token.
        """
        header = { "Authorization" : "Bearer %s" % (self.config["AUTH"]["access_token"], ) }
        return header

    def refresh_tokens(self):
        """
        Refresh any tokens nearing expiration.
        """
        now = int(time.time())
        any_changes = False
        for globus_service in ("AUTH", "TRANSFER", "GROUPS"):
            expires_in = int(self.config[globus_service]["expires_at_seconds"]) - now
            if DEBUG: print("Token %s expires in %s sec" % (globus_service, expires_in))
            if expires_in < 86400:
                any_changes = True
                refresh_token = self.config[globus_service]["refresh_token"]
                self._save_tokens_refresh(globus_service, globus_client.oauth2_refresh_token(refresh_token))
        if any_changes: self.config.write_new_config_file()

    def init(self):
        """
        Verify and refresh authorization.
        """
        self.refresh_tokens()
        self.scopes()
