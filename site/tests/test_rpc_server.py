import jaws_site.rpc_server


def test_rpc_start_server(mock_connection, mock_thread, mock_event):
    rpc_server = jaws_site.rpc_server.RpcServer()
    rpc_server.start_server()


def test_update_consumers(mock_connection, mock_thread, mock_event):
    rpc_server = jaws_site.rpc_server.RpcServer()
    rpc_server.increase_consumers_by(4)
    # our initial amount of consumers (defined as a default config) is 5
    # and so we add 4, thus total is 9 to test for.
    assert len(rpc_server.consumers) == 9
