"""Flask-SQLAlchemy db and models, used by Connexion/Flask server."""

import datetime
from flask_sqlalchemy import SQLAlchemy


db = SQLAlchemy()


class User(db.Model):
    """Registered user"""

    __tablename__ = "users"
    id = db.Column(db.String(32), primary_key=True)
    access_token = db.Column(db.String(256), nullable=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<User {self.id}>"
