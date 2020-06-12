"""
File that contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-central.ini"
    content = """[DB]
dialect = mysql+mysqlconnector
host = db.foo.com
port = 3306
user = jaws
password = passw0rd1
db = jaws

[GLOBUS]
client_id = ZZZZ

[SITE:LBNL]
host = rmq.foo.com
user = jaws
password = passw0rd2
vhost =
queue = jaws_rpc
globus_endpoint = XXXX
globus_basepath = "/global/scratch/jaws"
staging_subdir = "staging"
max_ram_gb = 1024

[SITE:NERSC]
host = rmq.bar.com
user = jaws
password = passw0rd3
vhost =
queue = jaws_rpc
globus_endpoint = YYYY
globus_basepath = "/"
staging_subdir = "/global/scratch/jaws/staging"
max_ram_gb = 2048
"""
    cfg.write_text(content)
    return cfg.as_posix()


@pytest.fixture()
def partial_config(tmp_path):
    cfg = tmp_path / "jaws-site.ini"
    content = """[DB]
dialect = mysql+mysqlconnector
host = db.foo.com
port = 3305
user = j4ws
password = p455w0rd1
db = jaws

[GLOBUS]
client_id = ZZZZ
"""
    cfg.write_text(content)
    return cfg.as_posix()


rpc_client_dict = {
    "host": "rmq.foo.com",
    "user": "jaws",
    "password": "passw0rd2",
    "vhost": "",
    "queue": "jaws_rpc",
    "globus_endpoint": "XXXX",
    "globus_basepath": '"/global/scratch/jaws"',
    "staging_subdir": "staging",
    "max_ram_gb": 1024,
}


@pytest.fixture()
def rpc_dict():
    return rpc_client_dict
