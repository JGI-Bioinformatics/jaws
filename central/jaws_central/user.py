from jaws_central import models_sa as models


class User:

    def __init__(self, session, uid: str):
        self.model = session.query(models.User).get(uid)

    @property
    def id(self):
        return self.model.id

    @property
    def is_admin(self):
        return self.model.is_admin
