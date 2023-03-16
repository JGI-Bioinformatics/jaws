import pytest
import jaws_site.rpc_operations
from jaws_site.runs import RunNotFoundError
from jaws_site import runs, transfers
from tests.conftest import (
    MockSession,
    MockRun,
    MockCromwell,
    MockCromwellException,
    MockTransferModel,
    MockTransfer,
)


def test_server_status(monkeypatch):
    monkeypatch.setattr(jaws_site.rpc_operations, "Cromwell", MockCromwell)
    mock_session = MockSession()
    p = {"user_id": "user", "run_id": 99}
    ret = jaws_site.rpc_operations.server_status(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "result": True}

    monkeypatch.setattr(jaws_site.rpc_operations, "Cromwell", MockCromwellException)
    mock_session = MockSession()
    p = {"user_id": "user", "run_id": 99}
    ret = jaws_site.rpc_operations.server_status(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "error": {"code": 500, "message": ""}}


def test_submit_run(monkeypatch):
    def mock_from_params(session, transfer_id):
        return MockRun()

    monkeypatch.setattr(runs.Run, "from_params", mock_from_params)

    mock_session = MockSession()
    p = {"user_id": "user", "run_id": 99}
    ret = jaws_site.rpc_operations.submit_run(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "result": "running"}

    def mock_from_params(session, transfer_id):
        raise Exception

    monkeypatch.setattr(runs.Run, "from_params", mock_from_params)

    ret = jaws_site.rpc_operations.submit_run(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "error": {"code": 500, "message": ""}}


@pytest.fixture()
def mock_run_from_id(monkeypatch):
    def mock_from_id(session, transfer_id):
        return MockRun()

    monkeypatch.setattr(runs.Run, "from_id", mock_from_id)


@pytest.fixture()
def mock_run_from_id_exception(monkeypatch):
    def mock_from_id(session, transfer_id):
        raise Exception

    monkeypatch.setattr(runs.Run, "from_id", mock_from_id)


@pytest.fixture()
def mock_run_from_id_runnotfound(monkeypatch):
    def mock_from_id(session, transfer_id):
        raise RunNotFoundError

    monkeypatch.setattr(runs.Run, "from_id", mock_from_id)


def test_run_all_functions(mock_run_from_id):
    mock_session = MockSession()
    p = {"user_id": "user", "run_id": 99}

    ret = jaws_site.rpc_operations.output_manifest(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "result": {"test": "success"}}

    ret = jaws_site.rpc_operations.cancel_run(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "result": {"test": "success"}}

    ret = jaws_site.rpc_operations.task_log(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "result": {"test": "success"}}


def test_run_all_functions_exception(mock_run_from_id_exception):
    mock_session = MockSession()
    p = {"user_id": "user", "run_id": 99}

    ret = jaws_site.rpc_operations.output_manifest(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "error": {"code": 500, "message": ""}}

    ret = jaws_site.rpc_operations.cancel_run(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "error": {"code": 500, "message": ""}}

    ret = jaws_site.rpc_operations.task_log(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "error": {"code": 500, "message": ""}}


def test_run_all_functions_runnotfound(mock_run_from_id_runnotfound):
    mock_session = MockSession()
    p = {"user_id": "user", "run_id": 99}

    ret = jaws_site.rpc_operations.output_manifest(p, mock_session)
    assert "error" in ret


@pytest.fixture()
def mock_transfer_from_id(monkeypatch):
    def mock_from_id(session, transfer_id):
        mock_session = MockSession()
        mock_data = MockTransferModel()
        return MockTransfer(mock_session, mock_data)

    monkeypatch.setattr(transfers.Transfer, "from_id", mock_from_id)


@pytest.fixture()
def mock_transfer_from_id_exception(monkeypatch):
    def mock_from_id(session, transfer_id):
        raise Exception

    monkeypatch.setattr(transfers.Transfer, "from_id", mock_from_id)


@pytest.fixture()
def mock_transfer_from_params(monkeypatch):
    def mock_from_params(session, transfer_id):
        mock_session = MockSession()
        mock_data = MockTransferModel()
        return MockTransfer(mock_session, mock_data)

    monkeypatch.setattr(transfers.Transfer, "from_params", mock_from_params)


@pytest.fixture()
def mock_transfer_from_params_exception(monkeypatch):
    def mock_from_params(session, transfer_id):
        raise Exception

    monkeypatch.setattr(transfers.Transfer, "from_params", mock_from_params)


def test_transfer_from_id_functions(mock_transfer_from_id):
    mock_session = MockSession()
    p = {
        "user_id": "user",
        "run_id": 99,
        "transfer_id": 123,
        "src_base_dir": "/test",
        "dest_base_dir": "/test",
    }

    ret = jaws_site.rpc_operations.transfer_status(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "result": {"status": "queued", "reason": None}}

    ret = jaws_site.rpc_operations.cancel_transfer(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "result": {"status": "canceled"}}


def test_transfer_from_id_functions_exception(mock_transfer_from_id_exception):
    mock_session = MockSession()
    p = {
        "user_id": "user",
        "run_id": 99,
        "transfer_id": 123,
        "src_base_dir": "/test",
        "dest_base_dir": "/test",
    }

    ret = jaws_site.rpc_operations.transfer_status(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "error": {"code": 500, "message": ""}}

    ret = jaws_site.rpc_operations.cancel_transfer(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "error": {"code": 500, "message": ""}}


def test_transfer_from_param_functions(mock_transfer_from_params):
    mock_session = MockSession()
    p = {
        "user_id": "user",
        "run_id": 99,
        "transfer_id": 123,
        "src_base_dir": "/test",
        "dest_base_dir": "/test",
    }

    ret = jaws_site.rpc_operations.submit_transfer(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "result": {"status": "queued"}}


def test_transfer_from_param_functions_exception(mock_transfer_from_params_exception):
    mock_session = MockSession()
    p = {
        "user_id": "user",
        "run_id": 99,
        "transfer_id": 123,
        "src_base_dir": "/test",
        "dest_base_dir": "/test",
    }

    ret = jaws_site.rpc_operations.submit_transfer(p, mock_session)
    assert ret == {"jsonrpc": "2.0", "error": {"code": 500, "message": ""}}
