import pytest
import os
from dataclasses import dataclass
from jaws_site.datatransfer_plugins import aws_transfer
# from jaws_site.datatransfer_protocol import DataTransferFactory, SiteTransfer


@pytest.fixture()
def mock_aws_transfer(monkeypatch):
    import boto3
    import threading

    @dataclass
    class Client:
        upload_file = "/my/upload_file"
        download_file = "/my/download_file"

    class MockupThread:
        @staticmethod
        def start(*args, **kwargs):
            pass

        @staticmethod
        def is_alive(*args, **kwargs):
            return data_obj.is_alive

    def mockup_client(*args, **kwargs):
        return Client

    def mockup_thread(*args, **kwargs):
        return MockupThread

    def mockup_get_file_path(self, src_path, dst_path):
        return src_path, dst_path

    monkeypatch.setattr(boto3, "client", mockup_client)
    monkeypatch.setattr(threading, "Thread", mockup_thread)
    monkeypatch.setattr(aws_transfer.DataTransfer, "_get_file_paths", mockup_get_file_path)

    class Data():
        def __init__(self):
            self.is_alive = True

    data_obj = Data()
    return data_obj


def test_add_transfer():
    result = aws_transfer.DataTransfer()._add_transfer()
    assert isinstance(result, str)


def test_submit_upload(mock_aws_transfer):
    metadata = {
        "label": "test",
    }
    manifest_files = [
        "/my/srcfile1\t/my/dstfile1\tF".encode('UTF-8'),
        "/my/srcfile2\t/my/dstfile2\tF".encode('UTF-8'),
        "/my/srcfile3\t/my/dstfile3\tF".encode('UTF-8'),
    ]
    result = aws_transfer.DataTransfer().submit_upload(metadata, manifest_files)
    assert result == 'True'


def test_submit_download(mock_aws_transfer):
    metadata = {
        "label": "test",
    }
    src_paths = [
        '/my/srcfile1',
        '/my/srcfile2',
        '/my/srcfile3',
    ]
    dst_paths = [
        '/my/dstfile1',
        '/my/dstfile2',
        '/my/dstfile3',
    ]
    result = aws_transfer.DataTransfer().submit_download(metadata, src_paths, dst_paths)
    assert result == 'True'


def test_cancel_transfer(mock_aws_transfer):
    result = aws_transfer.DataTransfer().cancel_transfer('123')
    assert result is None


def test_get_manifest_paths_filetype():
    # manifest list needs to be in bytestring as this is mimicking the return values from REST.
    manifest_files = [
        b"/my/srcfile1\t/my/dstfile1\tF",
        b"/my/srcfile2\t/my/dstfile2\tF",
        b"/my/srcfile3\t/my/dstfile3\tF",
    ]
    src_paths, dst_paths = aws_transfer.DataTransfer()._get_manifest_paths(manifest_files)
    assert src_paths == ['/my/srcfile1', '/my/srcfile2', '/my/srcfile3']
    assert dst_paths == ['/my/dstfile1', '/my/dstfile2', '/my/dstfile3']
    print(f"{src_paths=}")
    print(f"{dst_paths=}")


@pytest.fixture
def mock_src_paths(tmp_path):
    files = []
    f = tmp_path / "srcpath"
    f.mkdir()
    for n in range(1, 4):
        f = tmp_path / f"srcpath/srcfile{n}.txt"
        f.touch()
        files.append(f.as_posix())
    return files


def test_get_manifest_paths_dirtype(mock_src_paths):
    manifest_files = []
    for src_file in mock_src_paths:
        print(f"{src_file=}")
        src_dir = os.path.dirname(src_file)
        line = f"{src_dir}\t/my/dstfile\tD".encode('UTF-8')
        manifest_files.append(line)

    print(f"{manifest_files=}")
    src_paths, dst_paths = aws_transfer.DataTransfer()._get_manifest_paths(manifest_files)
    assert all(obs == exp for obs, exp in zip(src_paths, reversed(mock_src_paths)))
    assert dst_paths == ['/my/dstfile/srcfile3.txt', '/my/dstfile/srcfile2.txt', '/my/dstfile/srcfile1.txt']


def test_transfer_status(mock_aws_transfer):
    from jaws_site.datatransfer_protocol import SiteTransfer
    import threading

    # test transfer thread not defined, return failed status
    obj = aws_transfer.DataTransfer()
    mock_aws_transfer.is_alive = False
    result = obj.transfer_status('abcd')
    assert result == SiteTransfer.status.failed

    # test transfer thread defined, transferring still active, return tranferring status
    mock_aws_transfer.is_alive = True
    obj._transfer_threads['abcd'] = threading.Thread()
    result = obj.transfer_status('abcd')
    assert result == SiteTransfer.status.transferring

    # test transfer thread defined, transferring not active, return succeeded status
    mock_aws_transfer.is_alive = False
    result = obj.transfer_status('abcd')
    assert result == SiteTransfer.status.succeeded
