"""
File that contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-central.ini"
    content = """[JAWS]
name = jaws-dev
version = 2.0.1
docs_url = https://jaws-docs.readthedocs.io/en/latest/

[DB]
dialect = mysql+mysqlconnector
host = db.foo.com
port = 3306
user = jaws
password = passw0rd1
db = jaws

[GLOBUS]
client_id = ZZZZ

[SITE:LBNL]
host = rmq.jaws.gov
user = jaws
password = passw0rd2
vhost = jaws
queue = lbnl_rpc
port = 5672
globus_endpoint = XXXX
globus_basepath = "/global/scratch/jaws"
uploads_subdir = "uploads"
max_ram_gb = 1024

[SITE:NERSC]
host = rmq.jaws.gov
user = jaws
password = passw0rd2
vhost = jaws
queue = nersc_rpc
port = 5672
globus_endpoint = YYYY
globus_basepath = "/"
uploads_subdir = "/global/scratch/jaws/uploads"
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
    "host": "rmq.jaws.gov",
    "user": "jaws",
    "password": "passw0rd2",
    "vhost": "jaws",
    "queue": "lbnl_rpc",
    "globus_endpoint": "XXXX",
    "globus_basepath": '"/global/scratch/jaws"',
    "uploads_subdir": "uploads",
    "max_ram_gb": 1024,
}


@pytest.fixture()
def rpc_dict():
    return rpc_client_dict
