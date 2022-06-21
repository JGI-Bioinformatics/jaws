"""
File that contains all the mock classes and fixtures that will be used during
testing.
"""
import pytest
from datetime import datetime
import jaws_central


class MockSession:
    def __init__(self):
        return

    def commit(self):
        return

    def close(self):
        return

    def close_all(self):
        return

    def add(self, data):
        return

    def query(self, orm):
        return None


class MockDb:
    @property
    def session(self):
        return MockSession()


class MockRunModel:
    """Mock Run sqlalchemy orm model object with useable defaults."""

    def __init__(self, **kwargs):
        self.id = kwargs.get("id", "99")
        self.user_id = kwargs.get("user_id", "test_user")
        self.submission_id = kwargs.get("submission_id", "XXXX")
        self.caching = (
            False if "caching" in kwargs and kwargs["caching"] is False else True
        )
        self.input_site_id = kwargs.get("input_site_id", "CORI")
        self.compute_site_id = kwargs.get("compute_site_id", "JGI")
        self.cromwell_run_id = kwargs.get("cromwell_run_id", "myid")
        self.result = kwargs.get("result", "succeeded")
        self.status = kwargs.get("status", "running")
        self.submitted = kwargs.get("submitted", datetime.utcnow())
        self.updated = kwargs.get("updated", datetime.utcnow())
        self.upload_id = kwargs.get("upload_id", "12")
        self.download_id = kwargs.get("download_id", "13")
        self.wdl_file = kwargs.get("wdl_file", "/some/wdl")
        self.json_file = kwargs.get("json_file", "/some/json")
        self.tag = kwargs.get("tag", "some tag")
        self.manifest_json = kwargs.get("manifest_json", "{}")


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
inputs_dir = /global/scratch/jaws/jaws-dev/inputs
inputs_dir2 = ${SCRATCH_ROOT}/jaws/jaws-dev/inputs
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
inputs_dir = /global/cscratch/sd1/jaws/jaws-dev/inputs
max_ram_gb = 2048
[SITE:AWS]
host = rmq.jaws.gov
user = jaws
password = passw0rd4
vhost = jaws
queue = aws_rpc
message_ttl = 5
globus_endpoint =
globus_host_path =
inputs_dir = s3://jaws-site/jaws-dev/inputs
max_ram_gb = 512
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
    "inputs_dir": "/global/scratch/jaws/inputs",
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


class MockRunLogModel:

    def __init__(self, **kwargs):
        self.site_id = kwargs.get("site_id", "JGI")
        self.run_id = kwargs.get("run_id", 9527)
        self.status_from = kwargs.get("status_from", "queued")
        self.status_to = kwargs.get("status_to", "running")
        self.timestamp = kwargs.get("timestamp", datetime.utcnow())
        self.reason = kwargs.get("reason", None)


class MockTransferModel:
    """Mock Transfer sqlalchemy orm model object with useable defaults."""

    def __init__(self, **kwargs):
        self.id = kwargs.get("id", "99")
        self.status = kwargs.get("status", "created")
        self.src_site_id = kwargs.get("src_site_id", "NERSC")
        self.src_base_dir = kwargs.get("src_base_dir", "/jaws-test/inputs")
        self.dest_site_id = kwargs.get("dest_site_id", "JGI")
        self.dest_base_dir = kwargs.get("dest_base_dir", "/jaws-test/inputs")
        self.manifest_json = kwargs.get("manifest_json", "{}")
        self.globus_transfer_id = kwargs.get("globus_transfer_id", None)
        self.reason = kwargs.get("reason", None)


class MockTransfer:

    def __init__(self, session, data):
        pass

    @classmethod
    def from_id(cls, session, id):
        assert id is not None
        data = MockTransferModel(id=id)
        return cls(session, data)
