import jaws_site.rpc_server
import requests
import pytest
from tests.conftest import (
    MockMessage,
    mock_status_get,
    mock_metadata_get,
    mock_labels_get,
    mock_abort_get,
)


def test_rpc_start_server(
    amqp_config, rpc_config, mock_connection, mock_thread, mock_event
):
    rpc_server = jaws_site.rpc_server.RpcServer(amqp_config, rpc_config)
    rpc_server.start_server()


def test_update_consumers(
    amqp_config, rpc_config, mock_connection, mock_thread, mock_event
):
    rpc_server = jaws_site.rpc_server.RpcServer(amqp_config, rpc_config)
    rpc_server.update_consumers(0)
    rpc_server.update_consumers(10)
    assert 15 == len(rpc_server.consumers)


def test_dispatch_for_server_status(monkeypatch, mock_message):
    monkeypatch.setattr(requests, "get", mock_status_get)
    consumer = jaws_site.rpc_server.Consumer("test_queue")
    message = MockMessage("server_status", {}, "1")
    consumer(message)


def test_dispatch_for_fetching_metadata(monkeypatch, mock_message):
    monkeypatch.setattr(requests, "get", mock_metadata_get)
    consumer = jaws_site.rpc_server.Consumer("test_queue")
    message = MockMessage(
        "task_status", {"cromwell_id": "15774623-0f76-49ef-828c-3aa0ccd024f5"}, "1"
    )
    consumer(message)


def test_dispatch_for_fetching_labels(monkeypatch, mock_message):
    monkeypatch.setattr(requests, "get", mock_labels_get)
    consumer = jaws_site.rpc_server.Consumer("test_queue")
    message = MockMessage(
        "get_labels", {"cromwell_id": "15774623-0f76-49ef-828c-3aa0ccd024f5"}, "1"
    )
    consumer(message)


def test_dispatch_for_job_cancellation(monkeypatch, mock_message):
    monkeypatch.setattr(requests, "get", mock_abort_get)
    consumer = jaws_site.rpc_server.Consumer("test_queue")
    message = MockMessage(
        "cancel_run", {"cromwell_id": "15774623-0f76-49ef-828c-3aa0ccd024f5"}, "1"
    )
    consumer(message)


def test_raise_error_when_missing_cromwell_id(monkeypatch, mock_message):
    def mock_dispatch(method, params):
        """Return malformed JSON response"""
        return {"jsonrpc": "2.0", "error": []}

    monkeypatch.setattr(requests, "get", mock_metadata_get)
    monkeypatch.setattr(jaws_site.rpc_server, "dispatch", mock_dispatch)
    consumer = jaws_site.rpc_server.Consumer("test_queue")
    message = MockMessage("get_labels", {}, "1")
    with pytest.raises(jaws_site.rpc_server.InvalidJsonRpcResponse):
        consumer(message)
