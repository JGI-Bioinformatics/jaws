"""
File that contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-auth.ini"
    content = """[AUTH]
port = 3000

[DB]
dialect = mysql+mysqlconnector
host = db.foo.com
port = 3306
user = jaws
password = passw0rd1
db = jaws
"""
    cfg.write_text(content)
    return cfg.as_posix()


@pytest.fixture()
def partial_config(tmp_path):
    cfg = tmp_path / "jaws-auth.ini"
    content = """[DB]
host = db.foo.com
port = 3305
user = j4ws
password = p455w0rd1
db = jaws
"""
    cfg.write_text(content)
    return cfg.as_posix()


class MockDb:
    def __init__(self):
        return

    def session(self):
        return MockSession()


class MockSession:
    def __init__(self):
        return

    def commit(self):
        return

    def close_all(self):
        return

    def _query_user_id(self, *args, **kwargs):
        return


class MockUser:
    """Mock User object with useable defaults."""

    def __init__(self, **kwargs):
        self.id = kwargs.get("user_id", "Donald Duck")
        self.access_token = kwargs.get("access_token", "FFEEDDCCBBAA")
        self.is_admin = kwargs.get("is_admin", False)
