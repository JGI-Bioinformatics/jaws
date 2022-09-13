import jaws_site.rpc_operations
from jaws_site import runs
from tests.conftest import MockCromwellMetadata
from jaws_site import models
from datetime import datetime


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
    def __init__(self, session, data, raise_exception=False):
        self.session = session
        self.data = data
        self.raise_exception = raise_exception

    @classmethod
    def from_id(cls, session, id):
        assert id is not None
        data = MockTransferModel(id=id)
        return cls(session, data)

    def status(self):
        return self.data.status, self.data.reason

    def submit_transfer(self):
        if self.raise_exception:
            raise Exception
        pass

    def metadata(self):
        return MockCromwellMetadata("localhost/api/workflows/v1", "xxx-xxx")


def initRunModel(**kwargs):
    return models.Run(
        id=kwargs.get("id", "99"),
        user_id=kwargs.get("user_id", "test_user"),
        submission_id=kwargs.get("submission_id", "XXXX"),
        caching=(False if "caching" in kwargs and kwargs["caching"] is False else True),
        input_site_id=kwargs.get("input_site_id", "NERSC"),
        cromwell_run_id=kwargs.get("cromwell_run_id", "myid"),
        result=kwargs.get("result", "succeeded"),
        status=kwargs.get("status", "running"),
        submitted=kwargs.get("submitted", datetime.utcnow()),
        updated=kwargs.get("updated", datetime.utcnow()),
    )


def test_run_metadata(monkeypatch):
    def mock_from_id(session, transfer_id):
        mock_session = MockSession()
        return MockTransfer(mock_session, None)

    monkeypatch.setattr(runs.Run, "from_id", mock_from_id)

    mock_session = MockSession()
    p = {"user_id": "user", "run_id": 99}
    ret = jaws_site.rpc_operations.run_metadata(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "result": None}
