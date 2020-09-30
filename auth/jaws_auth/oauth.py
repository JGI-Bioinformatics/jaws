import logging
from flask import abort, request
from jaws_auth.api import get_user

logger = logging.getLogger(__package__)


def get_tokeninfo() -> dict:
    """
    OAuth2 REST endpoint: validate token in HTTP header and return user info dictionary.
    Abort on failure.
    """
    # Extract token from header.  Discard "Bearer" keyword by splitting
    try:
        _, access_token = request.headers["Authorization"].split()
    except Exception:
        abort(401, "Authentication failure; invalid header")

    # Get user record, if token exists
    try:
        user = get_user(access_token)
    except Exception as error:
        abort(500, f"Auth db unavailable: {error}; please try again later")

    # Abort if invalid token or return uid and scopes
    if user is None:
        logger.info(f"Invalid token {access_token}")
        abort(401, "Authentication failure")
    else:
        return user
