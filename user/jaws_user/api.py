"""
jaws-user classes
"""

import logging
import globus_sdk
from sqlalchemy.exc import SQLAlchemyError
from jaws_user import config
from jaws_user import db

logger = logging.getLogger(__package__)


class DatabaseError(Exception):
    pass


class UserNotFoundError(Exception):
    pass


class UserAlreadyExistsError(Exception):
    pass


class User:
    def __init__(self, user_id, params=None):
        """
        Get existing user or create new if params specified.
        """
        self.user_id = user_id
        if params:
            self.__add_user(params)
        else:
            self.__load_user(params)

    def __add_user(self):
        """
        Insert new user into the db.
        """
        logger.info(f"Add user {self.user_id}")
        session = db.Session()  # TODO parameterize
        try:
            user = session.query(db.User).get(self.user_id)
        except SQLAlchemyError as error:
            raise DatabaseError(f"{error}")
        if user:
            raise UserAlreadyExistsError()
        user = db.User(user_id=self.user_id)
        session.add(user)
        session.commit()
        session.close()

    def __load_user(self):
        """
        Retrieve an existing user's record from the db.
        """
        logger.debug(f"Get user {self.user_id}")
        session = db.Session()
        try:
            user = session.query(db.User).get(self.user_id)
        except SQLAlchemyError as error:
            session.rollback()
            session.close()
            raise DatabaseError(f"{error}")
        session.close()
        if user is None:
            raise UserNotFoundError()
        self.data = user

    def get_user(self):
        """
        Return a dict of the user info; it is a copy of the data, not the ORM object.
        This object cannot be manipulated by changing the dict.
        """
        user = self.data
        result = {
            "user_id": user.id,
            "email": user.email,
            "name": user.name,
            "is_admin": user.admin,
            "globus_id": user.globus_id,
            "auth_refresh_token": user.auth_refresh_token,
            "transfer_refresh_token": user.transfer_refresh_token,
        }
        return result

    def update_user(self, auth_refresh_token, transfer_refresh_token):
        """Save a user's Globus tokens.  Only the refresh tokens are required with the RefreshTokenAuthorizer
        and the refresh tokens never expire.  This is used when new users grant the jaws app globus access.

        :param auth_refresh_token: Auth Refresh Token
        :type auth_refresh_token: str
        :param auth_refresh_token: Transfer Refresh Token
        :type transfer_refresh_token: str
        """
        # get user info from Globus
        client = globus_sdk.NativeAppAuthClient(config.conf.get("GLOBUS", "client_id"))
        authorizer = globus_sdk.RefreshTokenAuthorizer(auth_refresh_token, client)
        auth_client = globus_sdk.AuthClient(authorizer=authorizer)
        # TODO add try block around Globus request
        user_info = auth_client.oauth2_userinfo()
        globus_id = user_info["sub"]
        name = user_info["name"]
        email = user_info["email"]

        # update row in db
        session = db.Session()  # TODO parameterize
        try:
            user = session.query(db.User).get(self.user_id)
        except SQLAlchemyError as error:
            raise DatabaseError(f"{error}")
        user.name = name
        user.email = email
        user.globus_id = globus_id
        user.auth_refresh_token = auth_refresh_token
        user.transfer_refresh_token = transfer_refresh_token
        try:
            session.commit()
        except Exception as error:
            session.rollback()
            session.close()
            raise DatabaseError(f"{error}")
        session.close()
