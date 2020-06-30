import jaws_rpc.rpc_server

test_config = {
    "host": "rmq.spaceforce.gov",
    "user": "blightyear",
    "password": "woody1",
    "vhost": "test",
    "queue": "test_rpc",
    "max_retries": "3",
    "num_threads": "5",
}


def test_method(params={}):
    return True


test_operations = {"test": test_method}


def test_rpc_start_server(mock_connection, mock_thread, mock_event):
    rpc_server = jaws_rpc.rpc_server.RpcServer(test_config, test_operations)
    rpc_server.start_server()


def test_update_consumers(mock_connection, mock_thread, mock_event):
    rpc_server = jaws_rpc.rpc_server.RpcServer(test_config, test_operations)
    rpc_server.increase_consumers_by(4)
    # our initial amount of consumers (defined as a default config) is 5
    # and so we add 4, thus total is 9 to test for.
    assert len(rpc_server.consumers) == 9


def test_dispatch(mock_connection, mock_thread, mock_event):
    rpc_server = jaws_rpc.rpc_server.RpcServer(test_config, test_operations)
    assert rpc_server.operations.get("test", {})
