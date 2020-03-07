import os
import sys
import requests
import globus_sdk
import configparser
import logging
from . import config

logger = None


class User:
    """
    User class
    """

    token_file = None
    globus_tokens = None
    _auth_client = None
    _transfer_client = None

    def __init__(self):
        """
        Load params from globus_tokens file
        """
        global logger
        logger = logging.getLogger(__package__)

        # LOAD TOKENS
        if "HOME" not in os.environ:
            sys.exit('Env var "HOME" not defined')
        self.token_file = os.path.join(
            os.environ["HOME"], f'.jaws.{config.conf.get("JAWS", "name")}.ini'
        )
        logger.info("Loading %s" % (self.token_file,))
        self.globus_tokens = configparser.ConfigParser()
        self.globus_tokens.read(self.token_file)

        # AUTHENTICATE AND INIT USER OBJECT
        access_token = self.globus_tokens["AUTH"]["access_token"]
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
        requested_scopes = (
            "openid profile email urn:globus:auth:scope:transfer.api.globus.org:all"
        )
        requested_scopes += (
            "urn:globus:auth:scope:groups.api.globus.org:view_my_groups_and_memberships"
        )
        GLOBUS_CLIENT_ID = config.conf.get("GLOBUS", "client_id")
        globus_client = globus_sdk.NativeAppAuthClient(GLOBUS_CLIENT_ID)
        globus_client.oauth2_start_flow(
            requested_scopes=requested_scopes, refresh_tokens=True
        )
        authorize_url = globus_client.oauth2_get_authorize_url()
        print(
            "Please go to this URL and log in\n(HINT: try command-click on the link):\n\n%s\n"
            % (authorize_url,)
        )
        auth_code = input("Then paste the authorization code here: ").strip()
        logger.info("Authenticating againt Globus")
        try:
            token_response = globus_client.oauth2_exchange_code_for_tokens(auth_code)
        except Exception:
            sys.exit("Authentication failed")

        # WRITE TOKENS TO USER CONFIG FILE
        auth_service_name = config.conf.get("GLOBUS", "auth_service_name")
        auth = token_response.by_resource_server[auth_service_name]
        self.globus_tokens["AUTH"]["access_token"] = auth["access_token"]
        self.globus_tokens["AUTH"]["refresh_token"] = auth["refresh_token"]
        self.globus_tokens["AUTH"]["expires_at_seconds"] = str(
            auth["expires_at_seconds"]
        )

        transfer_service_name = config.conf.get("GLOBUS", "transfer_service_name")
        transfer = token_response.by_resource_server[transfer_service_name]
        self.globus_tokens["TRANSFER"]["access_token"] = transfer["access_token"]
        self.globus_tokens["TRANSFER"]["refresh_token"] = transfer["refresh_token"]
        self.globus_tokens["TRANSFER"]["expires_at_seconds"] = str(
            transfer["expires_at_seconds"]
        )

        groups_service_name = config.conf.get("GLOBUS", "groups_service_name")
        groups = token_response.by_resource_server[groups_service_name]
        self.globus_tokens["GROUPS"]["access_token"] = groups["access_token"]
        self.globus_tokens["GROUPS"]["refresh_token"] = groups["refresh_token"]
        self.globus_tokens["GROUPS"]["expires_at_seconds"] = str(
            groups["expires_at_seconds"]
        )

        self.globus_tokens.write_new_token_file()

    def valid_token(self):
        """
        Validate JAWS token.  JAWS may return a new auth token, if the old one needed to be refreshed.
        """
        url = f'{config.conf.get("JAWS", "url")}/user'
        try:
            r = requests.get(url, headers=self.header())
        except Exception:
            sys.exit("Unable to communicate with JAWS server")
        if r.status_code != 200:
            sys.exit(r.text)

    #        result = r.json()
    #        if "new_access_token" in result:
    #            self.globus_tokens["AUTH"]["access_token"] = new_access_token
    #            self.globus_tokens.write_new_token_file()

    def header(self):
        """
        Return a dict to be used as HTTP OAuth2 header containing the current user's authentication token.
        """
        # PACK ALL GLOBUS TOKENS INTO A JAWS SUPER-TOKEN
        super_token = ":".join(
            (
                self.globus_tokens["AUTH"]["access_token"],
                self.globus_tokens["AUTH"]["refresh_token"],
                self.globus_tokens["AUTH"]["expires_at_seconds"],
                self.globus_tokens["TRANSFER"]["access_token"],
                self.globus_tokens["TRANSFER"]["refresh_token"],
                self.globus_tokens["TRANSFER"]["expires_at_seconds"],
                self.globus_tokens["GROUPS"]["access_token"],
                self.globus_tokens["GROUPS"]["refresh_token"],
                self.globus_tokens["GROUPS"]["expires_at_seconds"],
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
        access_token = self.globus_tokens["AUTH"]["access_token"]
        # NOTE: not using automatic refresh; see: refresh_tokens()
        # GLOBUS_CLIENT_ID = config.conf.get("GLOBUS", "client_id")
        # globus_client = globus_sdk.NativeAppAuthClient(GLOBUS_CLIENT_ID)
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
        access_token = self.globus_tokens["TRANSFER"]["access_token"]
        authorizer = globus_sdk.AccessTokenAuthorizer(access_token)
        self._transfer_client = globus_sdk.TransferClient(authorizer=authorizer)
        return self._transfer_client
