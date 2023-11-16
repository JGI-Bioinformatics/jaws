import json
import os
import os.path
import stat

import boto3
import botocore
import jaws_site
import pytest
import sqlalchemy
from deepdiff import DeepDiff
from jaws_site import transfers
from jaws_site.transfers import (FixPerms, Transfer, TransferDbError,
                                 TransferNotFoundError, TransferValueError,
                                 abs_to_rel_paths, check_fix_perms_queue,
                                 check_transfer_queue, get_abs_files,
                                 list_all_files_under_dir, mkdir,
                                 parallel_chmod)

from tests.conftest import (S3_BUCKET, MockSession, MockTransfer,
                            MockTransferModel, initTransferModel)


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
        self.TRANSFER_TYPE = "s3_download"

    def mock_s3_upload(self):
        self.TRANSFER_TYPE = "s3_upload"

    def mock_local_copy(self):
        self.TRANSFER_TYPE = "local_copy"

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
    transfer = Transfer(mock_session, mock_data)
    transfer.transfer_files()
    assert transfer.TRANSFER_TYPE == "s3_download"

    # if the dest path starts with "s3://" then upload to S3
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="/scratch/jaws/uploads/A/B/C",
        dest_base_dir="s3://jaws-site/uploads",
    )
    transfer = Transfer(mock_session, mock_data)
    transfer.transfer_files()
    assert transfer.TRANSFER_TYPE == "s3_upload"

    # if neither the src or dest path starts with "s3://" then rsync
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="/scratch/jaws/cromwell-executions/x",
        dest_base_dir="/scratch/jaws/downloads/x",
    )
    transfer = Transfer(mock_session, mock_data)
    transfer.transfer_files()
    assert transfer.TRANSFER_TYPE == "local_copy"


def test_local_copy(monkeypatch, mock_sqlalchemy_session, setup_files):
    src, dest = setup_files
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir=src,
        dest_base_dir=dest,
    )
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.local_copy()


def test_calculate_parallelism():
    max_threads = 10
    assert 1 == transfers.calculate_parallelism(10000)
    assert 5 == transfers.calculate_parallelism(50000)
    assert max_threads == transfers.calculate_parallelism(100000)
    assert max_threads == transfers.calculate_parallelism(10000000)
    assert max_threads == transfers.calculate_parallelism(335313)


def test_parallel_rsync_files_only(setup_files):
    src_base_dir, dest_base_dir = setup_files
    manifest = ["file99.txt"]
    jaws_site.transfers.parallel_rsync_files_only(manifest, src_base_dir, dest_base_dir)

    assert os.path.exists(os.path.join(dest_base_dir, "file99.txt"))
    for i in range(99):
        assert not os.path.exists(os.path.join(dest_base_dir, f"file{i}.txt"))


def test_handles_nonexistent_directory(mock_sqlalchemy_session):
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir="/src",
        dest_base_dir="/dst",
    )
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.transfer_files()
    assert transfer.data.status == "failed"


# @pytest.mark.skip(reason="this fails on mac but works on linux")
@pytest.mark.parametrize(
    "set_perms, expected_octal_perms", [("755", "0o755"), ("777", "0o777")]
)
def test_correctly_changes_permission(
    mock_sqlalchemy_session,
    monkeypatch,
    set_perms,
    config_file,
    expected_octal_perms,
    setup_files,
):
    monkeypatch.setenv("ENV_OVERRIDE_PREFIX", "ENV__")
    monkeypatch.setenv("ENV__SITE_file_permissions", set_perms)
    jaws_site.config.Configuration._destructor()
    conf = jaws_site.config.Configuration(
        config_file, env_prefix="ENV__"
    )  # recreate the config so we can override
    assert conf

    def get_permissions(path):
        return oct(stat.S_IMODE(os.stat(path).st_mode))

    src, dst = setup_files
    manifest = []
    for i in range(100):
        manifest.append(f"file{i}.txt")
    manifest_json = json.dumps(manifest)
    mock_data = MockTransferModel(
        status="queued",
        src_base_dir=src,
        dest_base_dir=dst,
        manifest_json=manifest_json,
    )
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.transfer_files()

    for i in range(100):
        assert (
            get_permissions(os.path.join(dst, f"file{i}.txt")) == expected_octal_perms
        )


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


def test_check_transfer_queue(mock_sqlalchemy_session, monkeypatch):
    mock_sqlalchemy_session.output(
        [
            {"id": 1},
            {"src_site_id": "NERSC"},
            {"src_base_dir": "/global/cscratch/jaws/jaws-dev/input"},
            {"dest_site_id": "AWS"},
            {"dest_base_dir": "s3://jaws-site/jaws-dev/inputs"},
        ]
    )
    check_transfer_queue(mock_sqlalchemy_session)
    assert mock_sqlalchemy_session.data.session["query"] is True

    monkeypatch.setattr(jaws_site.transfers, "Transfer", MockTransfer)

    # Test SQLAlchemyError
    mock_sqlalchemy_session.clear()
    mock_sqlalchemy_session.data.raise_exception_sqlalchemyerror = True
    check_transfer_queue(mock_sqlalchemy_session)


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


def test_s3_download(mock_sqlalchemy_session, monkeypatch):
    def mock_s3_download_files(self):
        self.TRANSFER_TYPE = "s3_download_files"

    def mock_s3_download_folder(self):
        self.TRANSFER_TYPE = "s3_download_folder"

    monkeypatch.setattr(transfers.Transfer, "s3_download_files", mock_s3_download_files)
    monkeypatch.setattr(
        transfers.Transfer, "s3_download_folder", mock_s3_download_folder
    )

    # TEST1: transfer files
    mock_data = initTransferModel(manifest_json='["file1"]')
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.s3_download()
    assert transfer.TRANSFER_TYPE == "s3_download_files"

    # TEST2: transfer folder
    mock_data = initTransferModel(manifest_json="[]")
    transfer = Transfer(mock_sqlalchemy_session, mock_data)
    transfer.s3_download()
    assert transfer.TRANSFER_TYPE == "s3_download_folder"


def test_s3_download_files(s3, mock_sqlalchemy_session, monkeypatch):
    def mock_manifest(self):
        return ["file1", "file2", "file3"]

    monkeypatch.setattr(transfers.Transfer, "manifest", mock_manifest)

    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)

    with pytest.raises(OSError):
        transfer.s3_download_files()


def test_s3_download_folder(s3, mock_sqlalchemy_session, monkeypatch):
    mock_data = initTransferModel()
    transfer = Transfer(mock_sqlalchemy_session, mock_data)

    with pytest.raises(botocore.exceptions.ParamValidationError):
        transfer.s3_download_folder()


def test_get_abs_files(setup_dir_tree):
    root = setup_dir_tree
    rel_paths = [
        "./file0.txt",
        "./a/file1.txt",
        "./a/b/file2.txt",
    ]
    expected = [
        f"{root}/file0.txt",
        f"{root}/a/file1.txt",
        f"{root}/a/b/file2.txt",
    ]
    actual = get_abs_files(root, rel_paths)
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False


def test_list_all_files_under_dir(setup_dir_tree):
    root = setup_dir_tree
    actual = list_all_files_under_dir(root)
    expected = [f"{root}/file0.txt", f"{root}/a/file1.txt", f"{root}/a/b/file2.txt"]
    assert actual == expected


def test_abs_to_rel_paths():
    root_dir = "/some/root"
    abs_paths = [
        "/some/root/file1",
        "/some/root/a/file2",
        "/some/root/a/b/file3",
    ]
    expected = ["file1", "a/file2", "a/b/file3"]
    actual = abs_to_rel_paths(root_dir, abs_paths)
    assert actual == expected


def test_parallel_chmod(setup_dir_tree):
    root = setup_dir_tree
    folders = [root, f"{root}/a", f"{root}/a/b"]
    files = [
        f"{root}/file0.txt",
        f"{root}/a/file1.txt",
        f"{root}/a/b/file2.txt",
    ]
    folder_mode_str = "0775"
    file_mode_str = "0664"
    folder_mode = int(folder_mode_str, base=8)
    file_mode = int(file_mode_str, base=8)
    parallelism = 2
    parallel_chmod(root, file_mode, folder_mode, parallelism)
    for folder in folders:
        actual = oct(os.stat(folder).st_mode)
        assert str(actual)[-4:] == folder_mode_str
    for file in files:
        actual = oct(os.stat(file).st_mode)
        assert str(actual)[-4:] == file_mode_str


def test_check_fix_perms_queue(mock_sqlalchemy_session, monkeypatch):
    def mock_fix_perms(self):
        assert self.data.base_dir is not None

    monkeypatch.setattr(FixPerms, "fix_perms", mock_fix_perms)
    mock_sqlalchemy_session.output(
        [
            {"id": 1, "base_dir": "/global/cscratch/jaws/jaws-dev/outputs/X"},
        ]
    )
    check_fix_perms_queue(mock_sqlalchemy_session)
    assert mock_sqlalchemy_session.data.session["query"] is True
