import pytest
import globus_sdk
from dataclasses import dataclass
from jaws_site.datatransfer_plugins import globus_transfer
from jaws_site.datatransfer_protocol import SiteTransfer, DataTransferError


@pytest.fixture()
def mock_globus_transfer(monkeypatch):
    class TransferClient:
        @staticmethod
        def get_task(*args, **kwargs):
            return {
                "status": data_obj.status
            }

        @staticmethod
        def submit_transfer(*args, **kwargs):
            result = {
                "task_id": "globus_task_id_1234"
            }
            return result

        @staticmethod
        def cancel_task(*args):
            return "cancelled"

    class TransferData:
        @staticmethod
        def add_item(*args, **kwargs):
            pass

    def mock_authorizer(*args, **kwargs):
        return

    def mock_transfer_client(*args, **kwargs):
        return TransferClient

    def mock_transfer_data(*args, **kwargs):
        return TransferData

    monkeypatch.setattr(globus_sdk, 'ClientCredentialsAuthorizer', mock_authorizer)
    monkeypatch.setattr(globus_sdk, 'TransferClient', mock_transfer_client)
    monkeypatch.setattr(globus_sdk, 'TransferData', mock_transfer_data)

    @dataclass
    class Data():
        status = "SUCCEEDED"  # other options: FAILED, INACTIVE

    data_obj = Data()
    return data_obj


def test_globus_submit_upload(mock_globus_transfer):
    metadata = {
        "label": "test",
        "host_paths": {
            "src": "/",
            "dest": "/"
        },
        "input_endpoint": 'abcd',
        "compute_endpoint": 'efgh',
        "run_id": 1234,
    }
    # manifest list needs to be in bytestring as this is mimicking the return values from REST.
    manifest_files = [
        b"/home/srcfile1\t/home/\tF",
        b"/home/srcfile2\t/home/\tF",
        b"/home/srcfile3\t/home/\tF",
    ]
    data_transfer = globus_transfer.DataTransfer()
    result = data_transfer.submit_upload(metadata, manifest_files)
    assert result == 'globus_task_id_1234'


def test_globus_submit_download(mock_globus_transfer):
    metadata = {
        "label": "test",
        "dest_endpoint": "abcd"
    }
    data_transfer = globus_transfer.DataTransfer()
    src_dir = '/global/scratch/jaws/test'
    dst_dir = '/home/destdir'
    result = data_transfer.submit_download(metadata, src_dir, dst_dir)
    assert result == 'globus_task_id_1234'


@pytest.mark.parametrize(
    'globus_status, expect_status',
    [
        ('FAILED', SiteTransfer.status.failed),
        ('INACTIVE', SiteTransfer.status.inactive),
        ('SUCCESS', SiteTransfer.status.succeeded),
    ]
)
def test_globus_transfer_status(globus_status, expect_status, mock_globus_transfer):
    mock_globus_transfer.status = globus_status
    data_transfer = globus_transfer.DataTransfer()
    status = data_transfer.transfer_status('123')
    assert status == expect_status


def test_globus_cancel_transfer(mock_globus_transfer):
    data_transfer = globus_transfer.DataTransfer()
    status = data_transfer.cancel_transfer('123')
    assert status == "cancelled"
