import logging
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from jaws_central import models


logger = logging.getLogger(__package__)


class UserError(Exception):
    # base class for all errors in this package
    pass


class UserDbError(UserError):
    pass


class UserNotFoundError(UserError):
    pass


class UserDataError(UserError):
    pass


class User:
    """Class representing a single User"""

    def __init__(self, session, data):
        self.session = session
        self.data = data

    @classmethod
    def from_params(cls, session, **kwargs):
        """Insert User record into rdb"""
        try:
            data = models.User(
                id=kwargs["id"],
                email=kwargs["email"],
                name=kwargs["name"],
                user_group=kwargs["user_group"],
                is_admin=True if kwargs["is_admin"] else False,
                is_dashboard=True if kwargs["is_dashboard"] else False,
                jaws_token=kwargs["jaws_token"]
            )
        except SQLAlchemyError as error:
            raise UserDataError(f"Error creating model for new User {kwargs['id']}: {error}")
        except Exception as error:
            raise UserDataError(f"Unknown error initializing User model: {error}")
        try:
            session.add(data)
            session.commit()
        except SQLAlchemyError as error:
            session.rollback()
            raise UserDbError(error)
        except Exception as error:
            session.rollback()
            raise UserDbError(error)
        else:
            return cls(session, data)

    @classmethod
    def from_id(cls, session, user_id):
        """Select user record from rdb"""
        try:
            data = session.query(models.User).get(user_id)
        except IntegrityError as error:
            logger.error(f"User {user_id} not found: {error}")
            raise UserNotFoundError(f"User {user_id} not found")
        except SQLAlchemyError as error:
            err_msg = f"Unable to select user, {user_id}: {error}"
            logger.error(err_msg)
            raise UserDbError(err_msg)
        except Exception as error:
            err_msg = f"Unexpected error; unable to select user, {user_id}: {error}"
            logger.error(err_msg)
            raise UserDbError(err_msg)
        else:
            return cls(session, data)
