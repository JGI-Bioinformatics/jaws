"""
JAWS User class provides authentication vs Globus and JAWS.
And provides the auth header method used by all requests vs JAWS Central.
"""

import os

# import grp
import sys
import requests
import globus_sdk

from jaws_client import config

# DEBUG = True if "JAWS_DEBUG" in os.environ else False
DEBUG = False
if "JAWS_DEBUG" in os.environ and os.environ["JAWS_DEBUG"]:
    DEBUG = True

# INIT GLOBUS PARAMS
JAWS_BRANCH = os.environ["JAWS_BRANCH"]
JAWS_URL = os.environ["JAWS_URL"]
GLOBUS_CLIENT_ID = os.environ["GLOBUS_CLIENT_ID"]
globus_client = globus_sdk.NativeAppAuthClient(GLOBUS_CLIENT_ID)
auth_service_name = "auth.globus.org"
transfer_service_name = "transfer.api.globus.org"
groups_service_name = (
    "04896e9e-b98e-437e-becd-8084b9e234a0"  # NOTE NOT groups.api.globus.org !?!
)


class User:
    """
    User class
    """

    config_file = None
    config = None
    _auth_client = None
    _transfer_client = None

    def __init__(self):
        """
        Load params from config file
        """
        # LOAD CONFIG
        if "HOME" not in os.environ:
            sys.exit('Env var "HOME" not defined')
        self.config_file = os.path.join(
            os.environ["HOME"], ".jaws.%s.ini" % (JAWS_BRANCH,)
        )
        if DEBUG:
            print("Loading %s" % (self.config_file,))
        self.config = config.Config(path=self.config_file, template="config.ini")

        # AUTHENTICATE AND INIT USER OBJECT
        access_token = self.config["AUTH"]["access_token"]
        if not access_token:
            print("Welcome, new user")
            if input("Do you have a Globus account? [y/n]: ").upper()[0] != "Y":
                sys.exit("You must create a Globus login first.")
            if (
                input("Have you been added to the JAWS group? [y/n]: ").upper()[0]
                != "Y"
            ):
                sys.exit(
                    "You must email an admin with your Globus username and request access to JAWS."
                )
            self.login()

    #        elif not self.valid_token():
    #            print("Your authentication has expired.")
    #            self.login()

    def login(self):
        """
        Authenticate using Globus OAuth2 and init object params.  Exit on failure.
        """
        # GLOBUS LOGIN
        requested_scopes = "openid profile email urn:globus:auth:scope:transfer.api.globus.org:all"
        requested_scopes += "urn:globus:auth:scope:groups.api.globus.org:view_my_groups_and_memberships"
        globus_client.oauth2_start_flow(
            requested_scopes=requested_scopes, refresh_tokens=True
        )
        authorize_url = globus_client.oauth2_get_authorize_url()
        print(
            "Please go to this URL and log in\n(HINT: try command-click on the link):\n\n%s\n"
            % (authorize_url,)
        )
        auth_code = input("Then paste the authorization code here: ").strip()
        if DEBUG:
            print("Authenticating againt Globus")
        try:
            token_response = globus_client.oauth2_exchange_code_for_tokens(auth_code)
        except Exception:
            sys.exit("Authentication failed")
        # if DEBUG: print(token_response)

        # WRITE TOKENS TO USER CONFIG FILE
        auth = token_response.by_resource_server[auth_service_name]
        self.config["AUTH"]["access_token"] = auth["access_token"]
        self.config["AUTH"]["refresh_token"] = auth["refresh_token"]
        self.config["AUTH"]["expires_at_seconds"] = str(auth["expires_at_seconds"])

        transfer = token_response.by_resource_server[transfer_service_name]
        self.config["TRANSFER"]["access_token"] = transfer["access_token"]
        self.config["TRANSFER"]["refresh_token"] = transfer["refresh_token"]
        self.config["TRANSFER"]["expires_at_seconds"] = str(
            transfer["expires_at_seconds"]
        )

        groups = token_response.by_resource_server[groups_service_name]
        self.config["GROUPS"]["access_token"] = groups["access_token"]
        self.config["GROUPS"]["refresh_token"] = groups["refresh_token"]
        self.config["GROUPS"]["expires_at_seconds"] = str(groups["expires_at_seconds"])

        if DEBUG:
            print("Saving token to file")
        self.config.write_new_config_file()

    def valid_token(self):
        """
        Validate JAWS token.  JAWS may return a new auth token, if the old one needed to be refreshed.
        """
        if DEBUG:
            print("Checking token with %s" % (self.header(),))
        url = "%s/user" % (JAWS_URL,)
        try:
            r = requests.get(url, headers=self.header())
        except Exception:
            sys.exit("Unable to communicate with JAWS server")
        if r.status_code != 200:
            sys.exit(r.text)
        if DEBUG:
            print("Welcome, back (valid token)")

    #        result = r.json()
    #        if "new_access_token" in result:
    #            self.config["AUTH"]["access_token"] = new_access_token
    #            self.config.write_new_config_file()

    def header(self):
        """
        Return a dict to be used as HTTP OAuth2 header containing the current user's authentication token.
        """
        # PACK ALL GLOBUS TOKENS INTO A JAWS SUPER-TOKEN
        super_token = ":".join(
            (
                self.config["AUTH"]["access_token"],
                self.config["AUTH"]["refresh_token"],
                self.config["AUTH"]["expires_at_seconds"],
                self.config["TRANSFER"]["access_token"],
                self.config["TRANSFER"]["refresh_token"],
                self.config["TRANSFER"]["expires_at_seconds"],
                self.config["GROUPS"]["access_token"],
                self.config["GROUPS"]["refresh_token"],
                self.config["GROUPS"]["expires_at_seconds"],
            )
        )
        header = {"Authorization": "Bearer %s" % (super_token,)}
        return header

    def auth_client(self):
        """
        Get Globus auth client object.
        """
        if self._auth_client is not None:
            return self._auth_client
        access_token = self.config["AUTH"]["access_token"]
        if DEBUG:
            print("User has auth token: %s" % (access_token,))
        # NOTE: not using automatic refresh; see: refresh_tokens()
        # authorizer = globus_sdk.RefreshTokenAuthorizer(refresh_token,
        #              globus_client, access_token=access_token, expires_at=expires_at_seconds)
        # authorizer = globus_sdk.RefreshTokenAuthorizer(refresh_token, globus_client)
        authorizer = globus_sdk.AccessTokenAuthorizer(access_token)
        self._auth_client = globus_sdk.AuthClient(authorizer=authorizer)
        return self._auth_client

    def transfer_client(self):
        """
        Get Globus transfer client object.
        """
        if self._transfer_client is not None:
            return self._transfer_client
        access_token = self.config["TRANSFER"]["access_token"]
        authorizer = globus_sdk.AccessTokenAuthorizer(access_token)
        self._transfer_client = globus_sdk.TransferClient(authorizer=authorizer)
        return self._transfer_client


#    def scopes(self):
#        """
#        Get JAWS scopes from Globus Groups and verify user is authorized to use JAWS.
#        """
#        if self._scopes is not None: return self._scopes
#        access_token = self.config["GROUPS"]["access_token"]
#        if not access_token: sys.exit("Globus Groups access token not defined")
#        headers={'Authorization': "Bearer %s" % (access_token,)}
#        r = requests.get(globus['groups_url'], headers=headers)
#        result = r.json()
#        scopes = []
#        for group in result:
#            group_name = group["name"]
#            # TODO use IDs not names
#            #group_id = group["id"]
#            #print("%s : %s" % (group_name, group_id)) # ECCE
#            if group_name.startswith("jaws_"): scopes.append(group_name)
#        self.scopes = scopes
#        if "jaws_users" not in scopes: sys.exit("You are not authorized to use JAWS")


#    def info(self):
#        """
#        Get user info from Globus.  The user's auxilary JAWS information is accessed separately (using SQLAlchemy).
#        """
#        if self._info is not None: return self._info
#        self._info = self.auth_client().oauth2_userinfo()
#        self._uid = self._info["sub"]
#        self._name = self._info["name"]
#        self._email = self._info["email"]
#        return self._info


#    def refresh_tokens(self):
#        """
#        Refresh any tokens nearing expiration.
#        """
#        epoch_time = int(time.time())
#        any_changes = False
#        for globus_service in ("AUTH", "TRANSFER", "GROUPS"):
#            expiry = int(self.config[globus_service]["expires_at_seconds"])
#            if expiry - epoch_time < 86400:
#                any_changes = True
#                refresh_token = self.config[globus_service]["refresh_token"]
#                self._save_tokens_refresh(globus_service, globus_client.oauth2_refresh_token(refresh_token))
#        if any_changes: self.config.write_new_config_file()
