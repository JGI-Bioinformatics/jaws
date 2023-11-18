import logging

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
    pass


test_operations = {"test": test_method}

dummy_sessionmaker = None


def test_rpc_start_server(mock_connection, mock_thread, mock_event):
    logger = logging.getLogger(__package__)
    rpc_server = jaws_rpc.rpc_server.RpcServer(
        test_config, logger, test_operations, dummy_sessionmaker
    )
    rpc_server.start_server()


def test_update_consumers(mock_connection, mock_thread, mock_event):
    logger = logging.getLogger(__package__)
    rpc_server = jaws_rpc.rpc_server.RpcServer(
        test_config, logger, test_operations, dummy_sessionmaker
    )
    rpc_server.increase_consumers_by(4)
    # our initial amount of consumers (defined as a default config) is 5
    # and so we add 4, thus total is 9 to test for.
    assert len(rpc_server.consumers) == 9


def test_dispatch(mock_connection, mock_thread, mock_event):
    logger = logging.getLogger(__package__)
    rpc_server = jaws_rpc.rpc_server.RpcServer(
        test_config, logger, test_operations, dummy_sessionmaker
    )
    assert rpc_server.operations.get("test", {})
