import pytest
import jaws_central.auth
import jaws_central.models_fsa

MOCK_JSON = {
    "ip": "1.2.3.4",
    "id": "abcd",
    "user": {
        "login": "JaneDoe",
        "email": "ijklmn@foo.gov",
        "first_name": "Jane",
        "middle_name": None,
        "last_name": "Doe",
        "email_address": "ijklmn@foo.gov"
    }
}


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

    def mock__get_json_from_sso(hash_code):
        return MOCK_JSON

    monkeypatch.setattr(
        jaws_central.auth, "_get_user_by_email", mock__get_user_by_email
    )

    monkeypatch.setattr(
        jaws_central.auth, "_get_json_from_sso", mock__get_json_from_sso
    )
    user = MockUser()
    hash_code = "abcde"
    result = jaws_central.auth.get_user_token(user, hash_code)

    assert result["jaws_token"] == "EEEEFFFFGGGG"
    assert result["sso_json"] == MOCK_JSON
    assert result["sso_json"]["user"]["first_name"] == "Jane"
