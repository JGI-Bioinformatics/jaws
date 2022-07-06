import pytest
import jaws_central.auth
import jaws_central.models
from tests.conftest import MockSession

MOCK_JSON = {
    "ip": "1.2.3.4",
    "id": "abcd",
    "user": {
        "login": "JaneDoe",
        "email": "ijklmn@foo.gov",
        "first_name": "Jane",
        "middle_name": None,
        "last_name": "Doe",
        "email_address": "ijklmn@foo.gov",
    },
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


@pytest.fixture()
def mock_database(monkeypatch):
    monkeypatch.setattr(jaws_central.models.db, "session", MockSession)


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
    # check if sending only email works and returns an empty SSO_JSON
    email = "ijklmn@foo.gov"
    result = jaws_central.auth.get_user_token(user, email)
    assert result["jaws_token"] == "EEEEFFFFGGGG"
    assert result["sso_json"] == {}

    sso_hash = "abcde"
    result = jaws_central.auth.get_user_token(user, sso_hash)
    assert result["jaws_token"] == "EEEEFFFFGGGG"
    assert result["sso_json"] == MOCK_JSON
    assert result["sso_json"]["user"]["first_name"] == "Jane"

    email = "abc.def@ghe@xyz"
    with pytest.raises(Exception) as e_info:
        result = jaws_central.auth.get_user_token(user, email)
    assert "Bad email address" in str(e_info.value)
