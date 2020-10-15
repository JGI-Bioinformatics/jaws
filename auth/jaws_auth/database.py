import flask_sqlalchemy


db = flask_sqlalchemy.SQLAlchemy()


class User(db.Model):
    """Registered user"""

    __tablename__ = "users"
    id = db.Column(db.String(32), primary_key=True)
    access_token = db.Column(db.String(256), nullable=False)
    is_admin = db.Column(db.Boolean(), nullable=False)
    # NOTE: defaults for Boolean columns fail due to type error with MySQL

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return f"<User {self.id}>"
