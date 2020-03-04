"""
SqlAlchemy ORM Models
"""

from sqlalchemy import (
    Column,
    String,
    Integer,
    Boolean,
    create_engine
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker

Base = declarative_base()


class User(Base):
    """
    Registered user record
    """

    __tablename__ = "users"
    id = Column(String(36), primary_key=True)
    name = Column(String(64), nullable=False)
    email = Column(String(64), nullable=False)
    is_admin = Column(Boolean, nullable=False, default=False)
    auth_access_token = Column(String(256))
    auth_refresh_token = Column(String(256))
    auth_expires_at_seconds = Column(Integer)
    transfer_access_token = Column(String(256))
    transfer_refresh_token = Column(String(256))
    transfer_expires_at_seconds = Column(Integer)
    groups_access_token = Column(String(256))
    groups_refresh_token = Column(String(256))
    groups_expires_at_seconds = Column(Integer)

    def update(
        self,
        name,
        email,
        auth_access_token,
        auth_refresh_token,
        auth_expires_at_seconds,
        transfer_access_token,
        transfer_refresh_token,
        transfer_expires_at_seconds,
        groups_access_token,
        groups_refresh_token,
        groups_expires_at_seconds,
    ):
        self.name = name
        self.email = email
        self.auth_access_token = auth_access_token
        self.auth_refresh_token = auth_refresh_token
        self.auth_expires_at_seconds = auth_expires_at_seconds
        self.transfer_access_token = transfer_access_token
        self.transfer_refresh_token = transfer_refresh_token
        self.transfer_expires_at_seconds = transfer_expires_at_seconds
        self.groups_access_token = groups_access_token
        self.groups_refresh_token = groups_refresh_token
        self.groups_expires_at_seconds = groups_expires_at_seconds


def init_db(uri):
    engine = create_engine(uri, convert_unicode=True)
    db_session = scoped_session(
        sessionmaker(autocommit=False, autoflush=False, bind=engine)
    )
    Base.query = db_session.query_property()
    Base.metadata.create_all(bind=engine)
    return db_session
