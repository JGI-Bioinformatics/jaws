"""Sqlalchemy ORM models for jaws-user database"""

from sqlalchemy import Column, String, Boolean
from jaws_user.database import Base


class User(Base):
    """Registered user"""

    __tablename__ = "users"
    id = Column(String(32), primary_key=True)
    email = Column(String(64), nullable=False)
    name = Column(String(64), nullable=True)
    is_admin = Column(Boolean, nullable=False, default=False)
    globus_id = Column(String(36), nullable=True)
    auth_refresh_token = Column(String(256), nullable=True)
    transfer_refresh_token = Column(String(256), nullable=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<User {self.id}>"
