"""
File that contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest
import jaws_central.config


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws-central.ini"
    content = """[JAWS]
name = jaws
version = 2.1
docs_url = https://jaws-docs.readthedocs.io/en/latest
[DB]
dialect = mysql+mysqlconnector
host = db.foo.com
host2 = ${JAWS_DB_HOST}
port = 3306
user = jaws
password = passw0rd1
password2 = ${JAWS_DB_PASSWORD}123
db = jaws
[RPC_SERVER]
host = rmq.jaws.gov
port = 5672
user = jaws
password = passw0rd2
vhost = jaws
queue = central_rpc
num_threads = 5
max_retries = 3
[HTTP]
auth_port = 3001
rest_port = 5001
[GLOBUS]
client_id = ZZZZ
client_secret = AAAAA
[SITE:JGI]
host = rmq.jaws.gov
user = jaws
password = passw0rd3
vhost = jaws
queue = lbnl_rpc
globus_endpoint = XXXX
globus_host_path = /global/scratch/jaws
uploads_dir = /global/scratch/jaws/jaws-dev/uploads
uploads_dir2 = ${SCRATCH_ROOT}/jaws/jaws-dev/uploads
max_ram_gb = 1024
[SITE:NERSC]
host = rmq.jaws.gov
user = jaws
password = passw0rd4
vhost = jaws
queue = nersc_rpc
message_ttl = 5
globus_endpoint = YYYY
globus_host_path = /
uploads_dir = /global/cscratch/sd1/jaws/jaws-dev/uploads
max_ram_gb = 2048
[ELASTIC_SEARCH]
host=jaws-vm-1.jgi.lbl.gob
port=9200
api_key=XXXX
[PERFORMANCE_METRICS]
index=perfmetrics
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

[RPC_SERVER]
user = jaws
password = pppaasss4
vhost = jaws_test

[GLOBUS]
client_id = AAAA
client_secret = BBBB
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
    "uploads_dir": "/global/scratch/jaws/uploads",
    "max_ram_gb": 1024,
}


@pytest.fixture()
def rpc_dict():
    return rpc_client_dict


@pytest.fixture()
def configuration(config_file):
    if jaws_central.config.conf is not None:
        jaws_central.config.Configuration._destructor()
    return jaws_central.config.Configuration(config_file)


@pytest.fixture()
def mock_data_transfer(monkeypatch):
    class MockDataTransfer:
        @staticmethod
        def submit_upload(metadata: str, manifest_file: list):
            return "123"

        @staticmethod
        def submit_download(metadata: dict, src_dir: str, dst_dir: str):
            return "456"

        @staticmethod
        def transfer_status(task_id: str):
            pass

        @staticmethod
        def cancel_transfer(task_id: str):
            return f"cancelled {task_id}"

    return MockDataTransfer
