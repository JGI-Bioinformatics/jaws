import logging
from flask import abort, request
from sqlalchemy.exc import SQLAlchemyError
from jaws_auth.models import db, User

logger = logging.getLogger(__package__)


def get_tokeninfo() -> dict:
    """
    OAuth2: validate token and return user info dictionary.
    Abort on login failure.
    """
    try:
        # Discard "Bearer" keyword by splitting
        _, access_token = request.headers["Authorization"].split()
    except Exception:
        abort(401, "Authentication failure; invalid header")

    try:
        user = (
            db.session.query(User).filter(User.access_token == access_token).one_or_none()
        )
    except SQLAlchemyError as e:
        abort(500, f"Db error: {e}")
    if user is None:
        logger.info(f"Authentication failure; got token {access_token}")
        abort(401, "Authentication failure")
    else:
        logger.debug(f"User {user.id} token OK")
        scopes = ["user"]
        return {"uid": user.id, "scope": scopes}
