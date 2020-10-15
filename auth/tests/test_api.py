import pytest
from deepdiff import DeepDiff
from jaws_auth import api


def test_get_user(monkeypatch):

    mock_table = {
        "AABBCCDDEEFF": {"id": "Mickey Mouse", "is_admin": True},
        "FFEEDDCCBBAA": {"id": "Donald Duck", "is_admin": False},
    }

    test_tokens = {
        "AABBCCDDEEFF": {"uid": "Mickey Mouse", "scopes": ["user", "admin"]},
        "FFEEDDCCBBAA": {"uid": "Donald Duck", "scopes": ["user"]},
    }

    def mock__select_from_db(self, session, token):
        if token not in mock_table:
            raise api.AuthenticationFailure()
        row = mock_table[token]
        return row["id"], row["is_admin"]

    monkeypatch.setattr(api.User, "_select_from_db", mock__select_from_db)

    # test with valid tokens
    for access_token in test_tokens:
        test_header = f"Bearer {access_token}"
        test_session = "NA"
        user = api.User(test_header, test_session)
        user_info = user.get_info()
        assert (
            bool(DeepDiff(user_info, test_tokens[access_token], ignore_order=True)) is False
        )

    # test with invalid token
    test_header = "Bearer THIS_TOKEN_DOES_NOT_EXIST"
    test_session = "NA"
    with pytest.raises(Exception):
        user = api.User(test_header, test_session)

    # test with invalid header
    test_header = "AABBCCDDEEFF"
    test_session = "NA"
    with pytest.raises(Exception):
        user = api.User(test_header, test_session)
