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
    def __init__(self, session=None):
        self.data = None
        if session:
            self.session = session
        else:
            self.session = db.session

    def create(self, user_id: str):
        """
        Insert new user into the db.
        """
        logger.info(f"Add user {self.user_id}")
        self.user_id = user_id

        # check if user exists
        try:
            user = self.session.query(db.User).get(user_id)
        except SQLAlchemyError as error:
            raise DatabaseError(f"{error}")
        if user:
            raise UserAlreadyExistsError()

        # insert
        try:
            self.data = db.User(user_id=self.user_id)
            self.session.add(user)
            self.session.commit()
        except Exception as error:
            self.session.rollback()
            raise DatabaseError(f"{error}")

    def load(self):
        """
        Retrieve an existing user's record from the db.
        """
        logger.debug(f"Get user {self.user_id}")
        try:
            self.data = self.session.query(db.User).get(self.user_id)
        except SQLAlchemyError as error:
            self.session.rollback()
            raise DatabaseError(f"{error}")
        if self.data is None:
            raise UserNotFoundError()

    def info(self):
        """
        Return a dict of the user info; it is a copy of the data, not the ORM object, so
        the object cannot be manipulated by changing the dict.
        """
        result = {
            "user_id": self.data.id,
            "email": self.data.email,
            "name": self.data.name,
            "is_admin": self.data.admin,
            "globus_id": self.data.globus_id,
            "auth_refresh_token": self.data.auth_refresh_token,
            "transfer_refresh_token": self.data.transfer_refresh_token,
        }
        return result

    def update(self, auth_refresh_token, transfer_refresh_token):
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
            raise DatabaseError(f"{error}")
