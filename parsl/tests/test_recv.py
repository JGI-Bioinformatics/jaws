import pytest
import jaws_parsl.recv
from unittest import mock

jaws_parsl.recv.rpc_params = {'user': 'jaws',
                              'password': 'p4s5w0rd',
                              'host': 'rpc.server.com',
                              'vhost': 'j4w5',
                              'port': 56789,
                              'queue': 'high-prio'
                              }


def test_start_file_logger(tmp_path):
    tmp_log = tmp_path / "test.txt"
    jaws_parsl.recv.start_file_logger(tmp_log)


def test_get_cromwell_run_id():
    msg = '/x/y/z/1-234-5/aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee/'
    expected_uuid = 'aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee'
    actual_uuid = jaws_parsl.recv.get_cromwell_run_id(msg)
    assert actual_uuid == expected_uuid


def test_get_cromwell_run_id_fail():
    msg = '/usr/bin/something'
    with pytest.raises(Exception):
        jaws_parsl.recv.get_cromwell_run_id(msg)


@mock.patch('jaws_parsl.recv.RpcClient', autospec=True)
def test_update_site(patch, monkeypatch):
    mock_rpc = mock.MagicMock()
    mock_rpc.__enter__.return_value = mock_rpc
    status_from = 'old'
    status_to = 'new'
    task_id = 1
    run_id = 'aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee'
    jaws_parsl.recv.update_site(status_from, status_to, task_id, run_id)


@mock.patch('jaws_parsl.recv.RpcClient', side_effect=Exception())
def test_update_site_fail(patch, monkeypatch):
    mock_rpc = mock.MagicMock()
    mock_rpc.__enter__.return_value = mock_rpc
    status_from = 'old'
    status_to = 'new'
    task_id = 1
    run_id = 'aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeeeeeee'
    with pytest.raises(Exception):
        jaws_parsl.recv.update_site(status_from, status_to, task_id, run_id)
