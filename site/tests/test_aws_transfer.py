import pytest
import os
from dataclasses import dataclass
from jaws_site.datatransfer_plugins import aws_transfer


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

    monkeypatch.setattr(boto3, "client", mockup_client)
    monkeypatch.setattr(threading, "Thread", mockup_thread)

    class Data():
        def __init__(self):
            self.is_alive = True

    data_obj = Data()
    return data_obj


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


def test_submit_download(monkeypatch, mock_aws_transfer):
    def mockup_get_file_paths(self, src_path, dst_path):
        return [src_path], [dst_path]

    monkeypatch.setattr(aws_transfer.DataTransfer, "_get_file_paths", mockup_get_file_paths)

    metadata = {
        "label": "test",
    }
    src_dir = '/global/scratch/jaws/test'
    dst_dir = '/home/destdir'
    result = aws_transfer.DataTransfer().submit_download(metadata, src_dir, dst_dir)
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


def test_get_manifest_paths_dirtype(mock_src_paths):
    manifest_files = []
    dst_dir = "/my/dstdir"
    for src_file in mock_src_paths:
        # src_file = <some_tmp_dir>/srcfile[1-3].txt created from mock_src_paths
        src_dir = os.path.dirname(src_file)
        line = f"{src_dir}\t{dst_dir}\tD".encode('UTF-8')
        manifest_files.append(line)

    src_paths, dst_paths = aws_transfer.DataTransfer()._get_manifest_paths(manifest_files)
    assert all(obs == exp for obs, exp in zip(src_paths, reversed(mock_src_paths)))
    assert dst_paths == [f"{dst_dir}/srcfile3.txt", f"{dst_dir}/srcfile2.txt", f"{dst_dir}/srcfile1.txt"]


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
