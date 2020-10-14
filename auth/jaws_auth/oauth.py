from flask import abort, request
from jaws_auth import api, db


def get_tokeninfo() -> dict:
    """
    OAuth2: validate token in HTTP header and return user info dictionary.
    Abort on failure.

    :return: user "uid" str and "scopes" list
    :rtype: dict
    """
    auth_header = request.headers["Authorization"]
    try:
        user = api.User(auth_header, db.session)
    except api.User.AuthInvalidHeader:
        abort(400, "Invalid header")
    except api.User.AuthDatabaseError:
        abort(500, "Auth db unavailable; please try again later")
    except api.User.AuthenticationError:
        abort(401, "Authentication failure")
    except Exception as error:
        abort(500, f"{error}")
    return user.get_info()
