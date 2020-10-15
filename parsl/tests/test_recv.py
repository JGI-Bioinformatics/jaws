import pytest
import jaws_parsl.recv
from unittest import mock
# from tests.__mocks__.jaws_rpc import rpc_client


# def test_update_site(monkeypatch):
#     mocked_rpc = Mock()
#     mocked_rpc.RPC_Client.return_value = rpc_client.RPC_Client()
#     monkeypatch.setattr('jaws_parsl.recv.RPC_Client', mocked_rpc)
#     jaws_parsl.recv.update_site('update', 1)

@mock.patch('jaws_parsl.recv.RPC_Client',autospec=True)
def test_update_site(patch):
    mock_am = mock.MagicMock()
    mock_am.__enter__.return_value = mock_am
    mock_am.request.return_value = {"a":"b"}
    jaws_parsl.recv.test_update_site('update', 1)


# @mock.patch('foobar.client_context_manager', autospec=True)
# def some_test(self, patch):
#     mock_client = mock.MagicMock(spec=DBClient)
#     mock_client.get_from_id = mock.Mock()
#     mock_client.foobar.side_effect = IndexError
#     mock_client.__enter__.return_value = mock_client
#     patch.return_value = mock_client