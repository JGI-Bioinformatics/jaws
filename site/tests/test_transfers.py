import pytest
import os
import os.path
import stat
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
    S3_BUCKET,
)
from subprocess import CalledProcessError
import sqlalchemy
import jaws_site
from jaws_site import transfers
import botocore
import boto3


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
    assert mock_session.needs_to_be_closed is False

    # a transfer that has already begun cannot be cancelled
    mock_data = MockTransferModel(status="transferring")
    transfer = Transfer(mock_session, mock_data)
    assert transfer.cancel() is False
    assert mock_session.needs_to_be_closed is False


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
    assert mock_session.needs_to_be_closed is False

    # if the dest path starts with "s3://" then upload to S3
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="/scratch/jaws/uploads/A/B/C",
        dest_base_dir="s3://jaws-site/uploads",
    )
    transfer = Transfer(mock_session, mock_data)
    transfer.transfer_files()
    assert transfer.S3_UPLOAD is True
    assert mock_session.needs_to_be_closed is False


def test_transfer_files2(mock_sqlalchemy_session, monkeypatch):
    def mock_rsync_folder(self):
        pass

    monkeypatch.setattr(Transfer, "rsync_folder", mock_rsync_folder)
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="/scratch/jaws-site/cromwell-executions/ex/AAAA",
        dest_base_dir="/scratch/jaws-site/uploads/AAAA",
    )
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.transfer_files()
    assert mock_sqlalchemy_session.data.session["commit"] is True
    assert transfer.data.status == "succeeded"


def test_calculate_parallelism():
    max_threads = 10
    assert 1 == transfers.calculate_parallelism(10000)
    assert 5 == transfers.calculate_parallelism(50000)
    assert max_threads == transfers.calculate_parallelism(100000)
    assert max_threads == transfers.calculate_parallelism(10000000)
    assert max_threads == transfers.calculate_parallelism(335313)


def test_use_find_subprocess_to_get_file_count(setup_files):
    src, _ = setup_files
    num_files = transfers.get_number_of_files(src)
    assert num_files == 1001  # Includes top-level directory


def test_parallel_rsync(mock_sqlalchemy_session, setup_files):
    src_base_dir, dest_base_dir = setup_files
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir=src_base_dir,
        dest_base_dir=dest_base_dir,
    )
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.transfer_files()

    for i in range(100):
        assert os.path.exists(os.path.join(dest_base_dir, f"file{i}.txt"))


def test_handles_nonexistent_directory(mock_sqlalchemy_session):
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="/src",
        dest_base_dir="/dst",
    )
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.transfer_files()
    assert transfer.data.status == "failed"


@pytest.mark.parametrize("set_perms, expected_octal_perms", [('775', "0o775"), ('777', "0o777")])
def test_correctly_changes_permission(mock_sqlalchemy_session, monkeypatch, set_perms,
                                             config_file, expected_octal_perms, setup_files):

    monkeypatch.setenv("ENV_OVERRIDE_PREFIX", "ENV__")
    monkeypatch.setenv("ENV__SITE_file_permissions", set_perms)
    jaws_site.config.Configuration._destructor()
    conf = jaws_site.config.Configuration(config_file, env_prefix="ENV__")  # recreate the config so we can override

    def get_permissions(path):
        return oct(stat.S_IMODE(os.stat(path).st_mode))

    src, dst = setup_files
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir=src,
        dest_base_dir=dst
    )
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.transfer_files()

    assert get_permissions(dst) == expected_octal_perms

    for i in range(100):
        assert get_permissions(os.path.join(dst, f"file{i}.txt")) == expected_octal_perms


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
        "close": True,
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
        "close": True,
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


def test_aws_s3_resource(s3, mock_sqlalchemy_session):
    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    s3_resource = transfer.aws_s3_resource()
    assert s3_resource is not None


def test_get_does_not_exist(s3, mock_sqlalchemy_session, monkeypatch):
    def mock_boto3_session(self):
        raise Exception

    monkeypatch.setattr(boto3, "Session", mock_boto3_session)

    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)

    with pytest.raises(Exception):
        transfer.aws_s3_resource()


class MockS3Session:
    def __init__(self):
        pass

    def resource(self, name):
        raise Exception


def test_aws_session_resource_exception(s3, mock_sqlalchemy_session, monkeypatch):

    monkeypatch.setattr(transfers.boto3, "Session", MockS3Session)

    with pytest.raises(Exception):
        mock_data = initTransferModel()
        transfer = Transfer(mock_sqlalchemy_session, mock_data)
        transfer.aws_s3_resource()


def test_aws_s3_client(s3, mock_sqlalchemy_session):
    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    cl = transfer.aws_s3_client()
    assert cl is not None

    response = s3.list_objects_v2(Bucket=S3_BUCKET)
    assert "Contents" not in response.keys()
    assert response["ResponseMetadata"]["HTTPStatusCode"] == 200


def test_s3_file_size(s3, mock_sqlalchemy_session, monkeypatch):
    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.s3_file_size(S3_BUCKET, "file_key")

    # Test list_objects_v2 ClientError exception
    class MockS3Client:
        def __init__(self):
            pass

        def list_objects_v2(self, Bucket=None, Prefix=None):
            raise botocore.exceptions.ClientError(
                error_response={"Error": {"Code": "code", "Message": "message"}},
                operation_name="operation_name",
            )

    def mock_aws_s3_client(self):
        return MockS3Client()

    monkeypatch.setattr(transfers.Transfer, "aws_s3_client", mock_aws_s3_client)

    with pytest.raises(botocore.exceptions.ClientError):
        transfer.s3_file_size(S3_BUCKET, "file_key")

    # Test list_objects_v2 ValueError exception
    class MockS3ClientParamValidationError:
        def __init__(self):
            pass

        def list_objects_v2(self, Bucket=None, Prefix=None):
            raise botocore.exceptions.ParamValidationError(report={})

    def mock_aws_s3_client(self):
        return MockS3ClientParamValidationError()

    monkeypatch.setattr(transfers.Transfer, "aws_s3_client", mock_aws_s3_client)

    with pytest.raises(ValueError):
        transfer.s3_file_size(S3_BUCKET, "file_key")

    # Test list_objects_v2  exception
    class MockS3ClientException:
        def __init__(self):
            pass

        def list_objects_v2(self, Bucket=None, Prefix=None):
            raise Exception

    def mock_aws_s3_client(self):
        return MockS3ClientException()

    monkeypatch.setattr(transfers.Transfer, "aws_s3_client", mock_aws_s3_client)

    with pytest.raises(Exception):
        transfer.s3_file_size(S3_BUCKET, "file_key")

    # Test if Contents, size=1
    class MockS3Client2:
        def __init__(self):
            pass

        def list_objects_v2(self, Bucket=None, Prefix=None):
            return {"Contents": [{"Size": 1}]}

    def mock_aws_s3_client(self):
        return MockS3Client2()

    monkeypatch.setattr(transfers.Transfer, "aws_s3_client", mock_aws_s3_client)

    transfer.s3_file_size(S3_BUCKET, "file_key")

    # Test if Contents, size>=1
    class MockS3Client3:
        def __init__(self):
            pass

        def list_objects_v2(self, Bucket=None, Prefix=None):
            return {
                "Contents": [
                    {"Size": 2, "Key": "file_key"},
                    {"Size": 3, "Key": "file_key"},
                ]
            }

    def mock_aws_s3_client(self):
        return MockS3Client3()

    monkeypatch.setattr(transfers.Transfer, "aws_s3_client", mock_aws_s3_client)

    transfer.s3_file_size(S3_BUCKET, "file_key")


def test_s3_upload(s3, mock_sqlalchemy_session, monkeypatch):
    def mock_manifest(self):
        return ["file1", "file2", "file3"]

    monkeypatch.setattr(transfers.Transfer, "manifest", mock_manifest)

    def mock_get_size(self):
        return 10

    monkeypatch.setattr(os.path, "getsize", mock_get_size)

    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)

    with pytest.raises(ValueError):
        transfer.s3_upload()


def test_s3_download(s3, mock_sqlalchemy_session, monkeypatch):
    def mock_manifest(self):
        return ["file1", "file2", "file3"]

    monkeypatch.setattr(transfers.Transfer, "manifest", mock_manifest)

    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)

    with pytest.raises(OSError):
        transfer.s3_download()


def test_s3_download_folder(s3, mock_sqlalchemy_session, monkeypatch):
    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)

    with pytest.raises(botocore.exceptions.ParamValidationError):
        transfer.s3_download_folder()
