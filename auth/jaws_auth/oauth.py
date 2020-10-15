from flask import abort, request
from jaws_auth.api import User, AuthInvalidHeader, AuthDatabaseError, AuthenticationError
from jaws_auth.database import db


def get_tokeninfo() -> dict:
    """
    OAuth2: validate token in HTTP header and return user info dictionary.
    Abort on failure.

    :return: user's "uid" str and "scopes" list
    :rtype: dict
    """
    auth_header = request.headers["Authorization"]
    try:
        user = User(auth_header, db.session)
    except AuthInvalidHeader:
        abort(400, "Invalid header")
    except AuthDatabaseError:
        abort(500, "Auth db unavailable; please try again later")
    except AuthenticationError:
        abort(401, "Authentication failure")
    except Exception as error:
        abort(500, f"{error}")
    return user.get_info()
