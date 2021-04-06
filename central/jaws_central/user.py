from jaws_central import db
from jaws_central.models_sa import User as UserModel

class User:

    def __init__(self, uid: str):
        self.model = db.session.query(UserModel).get(uid)

    @property
    def id(self):
        return self.model.id

    @property
    def is_admin(self):
        return self.model.is_admin
