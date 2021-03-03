import pytest
import jaws_central.auth
import jaws_central.models_fsa


class MockUser:
    @property
    def id(self):
        return "testuser"

    @property
    def jaws_token(self):
        return "EEEEFFFFGGGG"

    @property
    def email(self):
        return "ijklmn@foo.gov"

    @property
    def is_admin(self):
        return False

    @property
    def is_dashboard(self):
        return False


class MockQuery:
    @staticmethod
    def get(user):
        return MockUser()


class MockSession:
    @staticmethod
    def query(user):
        return MockQuery()


class MockDb:
    @property
    def session(self):
        return MockSession()


@pytest.fixture()
def mock_database(monkeypatch):
    monkeypatch.setattr(jaws_central.models_fsa.db, "session", MockSession)


def test_get_tokeninfo(monkeypatch):
    def mock__get_bearer_token():
        return {"Authorization": "Bearer AAAABBBBCCCC"}

    monkeypatch.setattr(jaws_central.auth, "_get_bearer_token", mock__get_bearer_token)

    def mock__get_user_by_token(access_token):
        user = MockUser()
        return user

    monkeypatch.setattr(
        jaws_central.auth, "_get_user_by_token", mock__get_user_by_token
    )

    result = jaws_central.auth.get_tokeninfo()

    assert result["uid"] == "testuser"


def test_get_user_token(monkeypatch):
    def mock__get_user_by_email(access_email):
        user = MockUser()
        return user

    monkeypatch.setattr(
        jaws_central.auth, "_get_user_by_email", mock__get_user_by_email
    )

    user = MockUser()
    email = "ijklmn@foo.gov"
    result = jaws_central.auth.get_user_token(user, email)

    assert result["jaws_token"] == "EEEEFFFFGGGG"
