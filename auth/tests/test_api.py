import pytest
from deepdiff import DeepDiff
from jaws_auth.api import get_user


def test_get_user(monkeypatch):

    mock_table = {
        "AABBCCDDEEFF": {"id": "Mickey Mouse", "is_admin": True},
        "FFEEDDCCBBAA": {"id": "Donald Duck", "is_admin": False},
    }

    test_tokens = {
        "AABBCCDDEEFF": {"uid": "Mickey Mouse", "scopes": ["user", "admin"]},
        "FFEEDDCCBBAA": {"uid": "Donald Duck", "scopes": ["user"]},
    }

    def mock___select_user_by_token(token):
        if token in mock_table:
            return mock_table[token]
        else:
            return None

    monkeypatch.setattr(auth, "__select_user_by_token", mock___select_user_by_token)

    for access_token in test_tokens:
        user = get_user(access_token)
        assert (
            bool(DeepDiff(user, test_tokens[access_token], ignore_order=True)) is False
        )
        user = get_user("this token doesn't exist")
        assert user is None
