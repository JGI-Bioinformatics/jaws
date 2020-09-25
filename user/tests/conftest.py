"""
File that contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-user.ini"
    content = """
[DB]
dialect = mysql+mysqlconnector
host = db.foo.com
port = 3306
user = jaws
password = passw0rd1
db = jaws

[RPC_SERVER]
user = jaws
password = pppaass4
vhost = jaws_test
"""
    cfg.write_text(content)
    return cfg.as_posix()


@pytest.fixture()
def partial_config(tmp_path):
    cfg = tmp_path / "jaws-user.ini"
    content = """[DB]
host = db.foo.com
port = 3305
user = j4ws
password = p455w0rd1
db = jaws

[RPC_SERVER]
user = jaws
password = pppaasss4
vhost = jaws_test
"""
    cfg.write_text(content)
    return cfg.as_posix()
