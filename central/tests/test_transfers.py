import pytest
from jaws_central import transfers
from tests.conftest import MockSession, MockTransferModel


def test_transfer_constructor():
    mock_session = MockSession()
    mock_data = MockTransferModel()
    transfer = transfers.Transfer(mock_session, mock_data)
    assert transfer


def test_status(monkeypatch):
    def mock_responsible_site_id_1(self):
        # None indicates no jaws-site is responsible (therefore check Globus)
        return None

    def mock_responsible_site_id_2(self):
        # string response is name of jaws-site to ask via RPC
        return "CORI"

    def mock_status_globus(self):
        return "transferring", None

    monkeypatch.setattr(transfers.Transfer, "status_globus", mock_status_globus)

    def mock_status_rpc(self, site_id):
        return "failed", "failure reason"

    monkeypatch.setattr(transfers.Transfer, "status_rpc", mock_status_rpc)

    def mock_update_status(self, new_status, reason):
        self.data.status = new_status
        self.data.reason = reason

    monkeypatch.setattr(transfers.Transfer, "update_status", mock_update_status)

    mock_session = MockSession()
    mock_data = MockTransferModel(status="queued")
    transfer = transfers.Transfer(mock_session, mock_data)

    monkeypatch.setattr(
        transfers.Transfer, "responsible_site_id", mock_responsible_site_id_1
    )
    actual, reason = transfer.status()
    assert actual == "transferring"

    monkeypatch.setattr(
        transfers.Transfer, "responsible_site_id", mock_responsible_site_id_2
    )
    actual, reason = transfer.status()
    assert actual == "failed"


def test_manifest():
    mock_session = MockSession()
    mock_data = MockTransferModel(manifest_json='{"foo": "bar"}')
    transfer = transfers.Transfer(mock_session, mock_data)
    manifest = transfer.manifest()
    assert type(manifest) is dict


def test_cancel():
    mock_session = MockSession()
    mock_data = MockTransferModel(status="created")
    transfer = transfers.Transfer(mock_session, mock_data)
    transfer.cancel()
    assert transfer.data.status == "cancelled"

    mock_data = MockTransferModel(status="completed")
    transfer = transfers.Transfer(mock_session, mock_data)
    transfer.cancel()
    assert transfer.data.status == "completed"


def test_responsible_site_id():
    mock_session = MockSession()

    # src and dest are the same
    test_src_site_id = "NERSC"
    test_dest_site_id = "NERSC"
    test_src_base_dir = "/global/cscratch/jaws/jaws-dev/inputs"
    test_dest_base_dir = "/global/cscratch/jaws/jaws-dev/inputs"
    mock_data = MockTransferModel(
        src_site_id=test_src_site_id,
        dest_site_id=test_dest_site_id,
        src_base_dir=test_src_base_dir,
        dest_base_dir=test_dest_base_dir,
    )
    transfer = transfers.Transfer(mock_session, mock_data)
    with pytest.raises(Exception):
        actual = transfer.responsible_site_id()

    # src and dest both have Globus endpoints
    test_src_site_id = "NERSC"
    test_dest_site_id = "JGI"
    test_src_base_dir = "/global/cscratch/jaws/jaws-dev/inputs"
    test_dest_base_dir = "/global/scratch/jaws/jaws-dev/inputs"
    mock_data = MockTransferModel(
        src_site_id=test_src_site_id,
        dest_site_id=test_dest_site_id,
        src_base_dir=test_src_base_dir,
        dest_base_dir=test_dest_base_dir,
    )
    transfer = transfers.Transfer(mock_session, mock_data)
    actual = transfer.responsible_site_id()
    assert actual is None

    # dest is AWS
    test_src_site_id = "NERSC"
    test_dest_site_id = "AWS"
    test_src_base_dir = "/global/cscratch/jaws/jaws-dev/inputs"
    test_dest_base_dir = "s3://jaws-site/jaws-dev/inputs"
    mock_data = MockTransferModel(
        src_site_id=test_src_site_id,
        dest_site_id=test_dest_site_id,
        src_base_dir=test_src_base_dir,
        dest_base_dir=test_dest_base_dir,
    )
    transfer = transfers.Transfer(mock_session, mock_data)
    actual = transfer.responsible_site_id()
    assert actual == test_src_site_id

    # src is AWS
    test_src_site_id = "AWS"
    test_dest_site_id = "NERSC"
    test_src_base_dir = "s3://jaws-site/jaws-dev/outputs"
    test_dest_base_dir = "/global/cscratch/jaws/jaws-dev/outputs"
    mock_data = MockTransferModel(
        src_site_id=test_src_site_id,
        dest_site_id=test_dest_site_id,
        src_base_dir=test_src_base_dir,
        dest_base_dir=test_dest_base_dir,
    )
    transfer = transfers.Transfer(mock_session, mock_data)
    actual = transfer.responsible_site_id()
    assert actual == test_dest_site_id
