import pytest
import os.path
from jaws_site.transfers import (
    Transfer,
    mkdir,
    TransferDbError,
    TransferValueError,
    TransferNotFoundError,
    check_queue,
)
from tests.conftest import (
    MockSession,
    MockTransferModel,
    initTransferModel,
    MockTransfer,
)
import sqlalchemy
import jaws_site


def test_mkdir(tmp_path):
    mkdir(tmp_path)
    assert os.path.exists(tmp_path)
    os.rmdir(tmp_path)


def test_constructor():
    mock_session = MockSession()
    mock_data = MockTransferModel()
    transfer = Transfer(mock_session, mock_data)
    assert transfer


def test_status():
    mock_session = MockSession()
    mock_data = MockTransferModel(status="queued")
    transfer = Transfer(mock_session, mock_data)
    assert transfer.status() == "queued"


def test_cancel(monkeypatch):
    def mock_update_status(self, new_status):
        assert type(new_status) is str
        assert new_status != self.data.status

    monkeypatch.setattr(Transfer, "update_status", mock_update_status)

    mock_session = MockSession()

    # a queued transfer may be cancelled
    mock_data = MockTransferModel(status="queued")
    transfer = Transfer(mock_session, mock_data)
    assert transfer.cancel() is True

    # a transfer that has already begun cannot be cancelled
    mock_data = MockTransferModel(status="transferring")
    transfer = Transfer(mock_session, mock_data)
    assert transfer.cancel() is False


def test_manifest():
    EXAMPLE_MANIFEST = ["file1", "file2", "file3"]
    mock_session = MockSession()
    mock_data = MockTransferModel(manifest=EXAMPLE_MANIFEST)
    transfer = Transfer(mock_session, mock_data)

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

    def mock_s3_download_folder(self):
        self.S3_DOWNLOAD_FOLDER = True

    def mock_s3_upload(self):
        self.S3_UPLOAD = True

    monkeypatch.setattr(Transfer, "s3_download", mock_s3_download)
    monkeypatch.setattr(Transfer, "s3_download_folder", mock_s3_download_folder)
    monkeypatch.setattr(Transfer, "s3_upload", mock_s3_upload)

    mock_session = MockSession()

    #    # if the src path starts with "s3://" then download from S3
    #    mock_data = MockTransferModel(
    #        status="queued",
    #        src_base_dir="s3://jaws-site/cromwell-executions/X/Y",
    #        dest_base_dir="/scratch/jaws/downloads",
    #    )
    #    transfer = Transfer(mock_session, mock_data)
    #    transfer.transfer_files()
    #    assert transfer.S3_DOWNLOAD is True

    # if the src path starts with "s3://" then download from S3
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="s3://jaws-site/cromwell-executions/AAAA",
        dest_base_dir="/scratch/jaws/downloads",
    )
    transfer = Transfer(mock_session, mock_data)
    transfer.transfer_files()
    assert transfer.S3_DOWNLOAD_FOLDER is True

    # if the dest path starts with "s3://" then upload to S3
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="/scratch/jaws/uploads/A/B/C",
        dest_base_dir="s3://jaws-site/uploads",
    )
    transfer = Transfer(mock_session, mock_data)
    transfer.transfer_files()
    assert transfer.S3_UPLOAD is True


def test_transfer_files2(mock_sqlalchemy_session, monkeypatch):
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="xxx://jaws-site/cromwell-executions/AAAA",
        dest_base_dir="xxx://jaws-site/uploads",
    )
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.transfer_files()
    # assert transfer.data.status == "failed"
    assert mock_sqlalchemy_session.data.session["commit"] is True
    assert transfer.data.status == "succeeded"

    # Test Exception
    # def mock_s3_download_folder(self):
    #     print('here')
    #     raise Exception
    #
    # monkeypatch.setattr(Transfer, "s3_download_folder", mock_s3_download_folder)
    #
    # mock_data = MockTransferModel(
    #     status="failed",
    #     src_base_dir="s3://jaws-site/cromwell-executions/AAAA",
    # )
    # transfer = Transfer(mock_sqlalchemy_session, mock_data)
    # transfer.transfer_files()
    # print(mock_sqlalchemy_session.data.session)


def test_s3_parse_path():
    mock_session = MockSession()
    mock_data = MockTransferModel()
    transfer = Transfer(mock_session, mock_data)

    test_path = "s3://jaws-site/dev/uploads/Pfam-A.hmm"
    actual_bucket, actual_path = transfer.s3_parse_path(test_path)
    assert actual_bucket == "jaws-site"
    assert actual_path == "dev/uploads/Pfam-A.hmm"


def test_from_params(mock_sqlalchemy_session):
    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    params = {
        "transfer_id": 99,
        "src_site_id": "NERSC",
        "dest_site_id": "AWS",
        "src_base_dir": "/global/cscratch/jaws/jaws-dev/inputs",
        "dest_base_dir": "s3://jaws-site/jaws-dev/inputs",
        "manifest": [],
    }
    transfer.from_params(mock_sqlalchemy_session, params)
    assert mock_sqlalchemy_session.data.session == {
        "add": True,
        "commit": True,
        "query": False,
        "close": False,
        "rollback": False,
    }
    mock_sqlalchemy_session.clear()

    params = {
        "transfer_id": 99,
        "src_site_id": "NERSC",
        "dest_site_id": "AWS",
        "src_base_dir": "/global/cscratch/jaws/jaws-dev/inputs",
        "dest_base_dir": "s3://jaws-site/jaws-dev/inputs",
        "manifest_json": "test",
    }
    transfer.from_params(mock_sqlalchemy_session, params)
    assert mock_sqlalchemy_session.data.session == {
        "add": True,
        "commit": True,
        "query": False,
        "close": False,
        "rollback": False,
    }
    mock_sqlalchemy_session.clear()

    # Test SQLAlchemyError
    mock_sqlalchemy_session.data.raise_exception_sqlalchemyerror = True
    with pytest.raises(TransferDbError):
        transfer.from_params(mock_sqlalchemy_session, params)
    assert mock_sqlalchemy_session.data.session["rollback"] is True
    mock_sqlalchemy_session.clear()

    params = {
        "transfer_id": "99",
        "src_base_dir": None,
        "dest_base_dir": None,
        "manifest_json": "test",
    }
    with pytest.raises(TransferValueError):
        transfer.from_params(mock_sqlalchemy_session, params)

    # Test wrong params
    params = {
        "src_base_dir": "/global/cscratch/jaws/jaws-dev/inputs",
        "dest_base_dir": "s3://jaws-site/jaws-dev/inputs",
        "manifest": [],
    }
    with pytest.raises(KeyError):
        transfer.from_params(mock_sqlalchemy_session, params)
    assert mock_sqlalchemy_session.data.session["rollback"] is True
    mock_sqlalchemy_session.clear()


def test_from_id(mock_sqlalchemy_session):
    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    mock_sqlalchemy_session.output(
        [
            {"src_site_id": "NERSC"},
            {"src_base_dir": "/global/cscratch/jaws/jaws-dev/input"},
            {"dest_site_id": "AWS"},
            {"dest_base_dir": "s3://jaws-site/jaws-dev/inputs"},
        ]
    )
    transfer.from_id(mock_sqlalchemy_session, 99)
    assert mock_sqlalchemy_session.data.session["query"] is True
    mock_sqlalchemy_session.clear()

    # Test Exception
    mock_sqlalchemy_session.data.raise_exception = True
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    with pytest.raises(Exception):
        transfer.from_id(mock_sqlalchemy_session, 99)
    mock_sqlalchemy_session.clear()

    # Test SQLAlchemyError
    mock_sqlalchemy_session.data.raise_exception_sqlalchemyerror = True
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    with pytest.raises(TransferDbError):
        transfer.from_id(mock_sqlalchemy_session, 99)
    mock_sqlalchemy_session.clear()

    # None result from one_or_none()
    mock_sqlalchemy_session.output([])
    with pytest.raises(TransferNotFoundError):
        transfer.from_id(mock_sqlalchemy_session, 99)


def test_reason(mock_sqlalchemy_session):
    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    ret = transfer.reason()
    assert ret == "reason"


def test_update_status(mock_sqlalchemy_session):
    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.update_status("queued", "")
    assert mock_sqlalchemy_session.data.session["commit"] is True
    mock_sqlalchemy_session.clear()

    # Test TransferDbError
    mock_sqlalchemy_session.data.raise_exception_sqlalchemyerror = True
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    with pytest.raises(sqlalchemy.exc.SQLAlchemyError):
        transfer.update_status("queued", "")
    assert mock_sqlalchemy_session.data.session["rollback"] is True
    mock_sqlalchemy_session.clear()


def test_check_queue(mock_sqlalchemy_session, monkeypatch):
    mock_sqlalchemy_session.output(
        [
            {"id": 1},
            {"src_site_id": "NERSC"},
            {"src_base_dir": "/global/cscratch/jaws/jaws-dev/input"},
            {"dest_site_id": "AWS"},
            {"dest_base_dir": "s3://jaws-site/jaws-dev/inputs"},
        ]
    )
    check_queue(mock_sqlalchemy_session)
    assert mock_sqlalchemy_session.data.session["query"] is True

    monkeypatch.setattr(jaws_site.transfers, "Transfer", MockTransfer)

    # Test SQLAlchemyError
    mock_sqlalchemy_session.clear()
    mock_sqlalchemy_session.data.raise_exception_sqlalchemyerror = True
    check_queue(mock_sqlalchemy_session)
