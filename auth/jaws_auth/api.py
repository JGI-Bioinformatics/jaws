import logging
from sqlalchemy.exc import SQLAlchemyError
from jaws_auth import db

logger = logging.getLogger(__package__)


class AuthInvalidHeader(Exception):
    pass


class AuthDatabaseError(Exception):
    pass


class AuthenticationFailure(Exception):
    pass


class User:
    def __init__(self, auth_header: str, session: db.Session = None):
        """
        Get user from db, if token exists.  Abort otherwise.

        :param token: HTTP authorization header
        :type token: str
        :return: user record or None if user not found
        :rtype: dict
        """
        if not session:
            try:
                session = db.session()
            except SQLAlchemyError as error:
                logger.exception(f"Unable to get session: {error}")
                raise AuthDatabaseError()

        # Extract token from header.  Discard "Bearer" keyword by splitting
        try:
            _, access_token = auth_header.split()
        except Exception:
            raise AuthInvalidHeader()

        # query db
        try:
            row = (
                self.session.query(db.User)
                .filter(User.access_token == access_token)
                .one_or_none()
            )
        except SQLAlchemyError as error:
            logger.exception(f"Select failed: {error}")
            raise AuthDatabaseError()

        # raise if invalid token
        if not row:
            logger.warn(f"Invalid token {access_token}")
            raise AuthenticationFailure()

        self.data = row

    def get_info(self):
        """
        Get info about current user.

        :return: uid (str) and scopes (list):
        :rtype: dict
        """
        user_info = {
            "uid": self.data.id,
            "scopes": ["user"]
        }
        if self.data.is_admin:
            user_info["scopes"].append("admin")
        return user_info
