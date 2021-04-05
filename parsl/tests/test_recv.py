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


@mock.patch('jaws_parsl.recv.RpcClient', autospec=True)
def test_update_site(patch, monkeypatch):
    mock_rpc = mock.MagicMock()
    mock_rpc.__enter__.return_value = mock_rpc
    jaws_parsl.recv.update_site('update', 1)


@mock.patch('jaws_parsl.recv.RpcClient', side_effect=Exception())
def test_update_site_fail(patch, monkeypatch):
    mock_rpc = mock.MagicMock()
    mock_rpc.__enter__.return_value = mock_rpc
    with pytest.raises(Exception):
        jaws_parsl.recv.update_site('update', 1)
