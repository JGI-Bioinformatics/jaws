import pytest
import jaws_parsl.recv
from unittest import mock

@mock.patch('jaws_parsl.recv.RpcClient',autospec=True)
def test_update_site(patch, monkeypatch):
    mock_rpc = mock.MagicMock()
    mock_rpc.__enter__.return_value = mock_rpc
    mock_rpc.request.return_value = {"jsonrpc": "2.0", "result": {"foo": "bar"}}
    mock_log = mock.Mock()
    monkeypatch.setattr('jaws_parsl.recv.logger', mock_log)
    jaws_parsl.recv.update_site('update', 1)


@mock.patch('jaws_parsl.recv.RpcClient', side_effect=Exception())
def test_update_site_fail(patch, monkeypatch):
    mock_rpc = mock.MagicMock()
    mock_rpc.__enter__.return_value = mock_rpc
    mock_log = mock.Mock()
    monkeypatch.setattr('jaws_parsl.recv.logger', mock_log)
    with pytest.raises(Exception):
        jaws_parsl.recv.update_site('update', 1)
