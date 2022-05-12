import shutil
import logging
from jaws_site import transfers
from jaws_site.transfers import Transfer
from tests.conftest import MockSession, MockTransferModel


logger = logging.getLogger(__package__)


def test_constructor():
    mock_session = MockSession()
    mock_data = MockTransferModel()
    transfer = Transfer(mock_session, logger, mock_data)
    assert transfer


def test_status():
    mock_session = MockSession()
    mock_data = MockTransferModel(status="queued")
    transfer = Transfer(mock_session, logger, mock_data)
    assert transfer.status() == "queued"


def test_cancel(monkeypatch):
    def mock_update_status(self, new_status):
        assert type(new_status) is str
        assert new_status != self.data.status

    monkeypatch.setattr(Transfer, "update_status", mock_update_status)

    mock_session = MockSession()

    # a queued transfer may be cancelled
    mock_data = MockTransferModel(status="queued")
    transfer = Transfer(mock_session, logger, mock_data)
    assert transfer.cancel() is True

    # a transfer that has already begun cannot be cancelled
    mock_data = MockTransferModel(status="transferring")
    transfer = Transfer(mock_session, logger, mock_data)
    assert transfer.cancel() is False


def test_manifest():

    EXAMPLE_MANIFEST = ["file1", "file2", "file3"]
    mock_session = MockSession()
    mock_data = MockTransferModel(manifest=EXAMPLE_MANIFEST)
    transfer = Transfer(mock_session, logger, mock_data)

    assert type(transfer.data.manifest_json) == str
    manifest = transfer.manifest()
    assert type(manifest) == list
    assert type(manifest[0]) == str


def test_transfer_files(monkeypatch):
    def mock_update_status(self, new_status):
        assert type(new_status) is str
        assert new_status != self.data.status

    monkeypatch.setattr(Transfer, "update_status", mock_update_status)

    # these 3 mocks set a variable in the object so we can see which method was called

    def mock_s3_download(self):
        self.S3_DOWNLOAD = True

    def mock_s3_upload(self):
        self.S3_UPLOAD = True

    def mock_local_copy(self):
        self.LOCAL_COPY = True

    monkeypatch.setattr(Transfer, "s3_download", mock_s3_download)
    monkeypatch.setattr(Transfer, "s3_upload", mock_s3_upload)
    monkeypatch.setattr(Transfer, "local_copy", mock_local_copy)

    mock_session = MockSession()

    # if the src path starts with "s3://" then download from S3
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="s3://jaws-site/cromwell-executions/X/Y",
        dest_base_dir="/scratch/jaws/downloads",
    )
    transfer = Transfer(mock_session, logger, mock_data)
    transfer.transfer_files()
    assert transfer.S3_DOWNLOAD is True

    # if the dest path starts with "s3://" then upload to S3
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="/scratch/jaws/uploads/A/B/C",
        dest_base_dir="s3://jaws-site/uploads",
    )
    transfer = Transfer(mock_session, logger, mock_data)
    transfer.transfer_files()
    assert transfer.S3_UPLOAD is True

    # otherwise local copy
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="/scratch/jaws/cromwell-executions/X/Y",
        dest_base_dir="/scratch/jaws/downloads",
    )
    transfer = Transfer(mock_session, logger, mock_data)
    transfer.transfer_files()
    assert transfer.LOCAL_COPY is True


def test_local_copy(monkeypatch):
    def mock_update_status(self, new_status):
        assert type(new_status) is str
        assert new_status != self.data.status

    def mock_mkdir(path):
        assert path.startswith("/")

    def mock_copyfile(self, src_path, dest_path):
        assert src_path != dest_path
        assert src_path.startswith("/scratch/jaws-/cromwell-executions")
        assert dest_path.startswith("/scratch/jaws/downloads")

    monkeypatch.setattr(Transfer, "update_status", mock_update_status)
    monkeypatch.setattr(transfers, "mkdir", mock_mkdir)
    monkeypatch.setattr(shutil, "copyfile", mock_copyfile)

    mock_session = MockSession()

    # if the src path starts with "s3://" then download from S3
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="/scratch/jaws/cromwell-executions",
        dest_base_dir="/scratch/jaws/downloads",
        manifest_json='{"task1.file1": "file1.txt"}',
    )
    transfer = Transfer(mock_session, logger, mock_data)
    transfer.local_copy()


def test_s3_parse_path():
    mock_session = MockSession()
    mock_data = MockTransferModel()
    transfer = Transfer(mock_session, logger, mock_data)

    test_path = "s3://jaws-site/dev/uploads/Pfam-A.hmm"
    actual_bucket, actual_path = transfer.s3_parse_path(test_path)
    assert actual_bucket == "jaws-site"
    assert actual_path == "dev/uploads/Pfam-A.hmm"
