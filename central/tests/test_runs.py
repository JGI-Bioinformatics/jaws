from jaws_central import runs, transfers
from tests.conftest import MockSession, MockRunModel, MockTransferModel


def test_run_constructor():
    mock_session = MockSession()
    mock_data = MockRunModel()
    run = runs.Run(mock_session, mock_data)
    assert run


def test_run_operations():
    mock_session = MockSession()
    mock_data = MockRunModel()
    run = runs.Run(mock_session, mock_data)

    required_operations = [
        "created",
        "upload queued",
        "uploading",
        "upload complete",
        "finished",
        "download queued",
        "downloading",
    ]
    for required_operation in required_operations:
        assert required_operation in run.operations
        assert callable(run.operations[required_operation])


def test_get_upload(monkeypatch):
    def mock_get_transfer(self, transfer_id):
        session = MockSession()
        data = MockTransferModel(id=transfer_id)
        transfer = transfers.Transfer(session, data)
        return transfer

    monkeypatch.setattr(runs.Run, "_get_transfer", mock_get_transfer)

    session = MockSession()

    data = MockRunModel(upload_id=None)
    run = runs.Run(session, data)
    transfer = run.get_upload()
    assert transfer is None

    test_transfer_id = 999
    data = MockRunModel(upload_id=test_transfer_id)
    run = runs.Run(session, data)
    transfer = run.get_upload()
    assert transfer.data.id == test_transfer_id


def test_get_download(monkeypatch):
    def mock_get_transfer(self, transfer_id):
        session = MockSession()
        data = MockTransferModel(id=transfer_id)
        transfer = transfers.Transfer(session, data)
        return transfer

    monkeypatch.setattr(runs.Run, "_get_transfer", mock_get_transfer)

    session = MockSession()

    data = MockRunModel(download_id=None)
    run = runs.Run(session, data)
    transfer = run.get_download()
    assert transfer is None

    test_transfer_id = 999
    data = MockRunModel(download_id=test_transfer_id)
    run = runs.Run(session, data)
    transfer = run.get_download()
    assert transfer.data.id == test_transfer_id


def test_check_if_upload_complete(monkeypatch):
    def mock_get_transfer(self, transfer_id):
        session = MockSession()
        data = MockTransferModel(id=transfer_id, status="succeeded")
        transfer = transfers.Transfer(session, data)
        return transfer

    def mock_update_status(self, new_status):
        self.data.status = new_status

    monkeypatch.setattr(runs.Run, "_get_transfer", mock_get_transfer)
    monkeypatch.setattr(runs.Run, "update_status", mock_update_status)

    session = MockSession()
    test_transfer_id = 999
    data = MockRunModel(upload_id=test_transfer_id, status="uploading")
    run = runs.Run(session, data)
    assert run.data.status == "uploading"
    run.check_if_upload_complete()
    assert run.data.status == "upload complete"


def test_check_if_download_complete(monkeypatch):
    def mock_get_transfer(self, transfer_id):
        session = MockSession()
        data = MockTransferModel(id=transfer_id, status="failed", reason='failure reason')
        transfer = transfers.Transfer(session, data)
        return transfer

    def mock_update_status(self, new_status, reason):
        self.data.status = new_status
        self.data.reason = reason

    monkeypatch.setattr(runs.Run, "_get_transfer", mock_get_transfer)
    monkeypatch.setattr(runs.Run, "update_status", mock_update_status)

    session = MockSession()
    test_transfer_id = 999
    data = MockRunModel(download_id=test_transfer_id, status="downloading")
    run = runs.Run(session, data)
    assert run.data.status == "downloading"
    run.check_if_download_complete()
    assert run.data.status == "download failed"


# def test_cancel_run(monkeypatch):
#    """Test the cancel run functioning."""
#
#    def get_run_cancelled(user_id, run_id):
#        run = MockGetInactiveUploadRun()
#        run.status = "cancelled"
#        return run
#
#    def get_run_regular(user_id, run_id):
#        run = MockGetInactiveUploadRun()
#        run.status = "running"
#        return run
#
#    def mock_cancel_run(user, run):
#        run.status = "cancelled"
#        run.result = "cancelled"
#        return run
#
#    def mock_cancel_transfer(transfer_task_id):
#        pass
#
#    def mock_rpc_call_cancel(user_id, run_id, method, params={}):
#        assert isinstance(user_id, str)
#        assert isinstance(run_id, int)
#        assert method == "cancel_run"
#
#    """Check if an exception is raised in case run is already cancelled"""
#    monkeypatch.setattr(jaws_central.rest, "rpc_call", mock_rpc_call_cancel)
#    monkeypatch.setattr(jaws_central.rest, "_get_run", get_run_cancelled)
#    monkeypatch.setattr(jaws_central.rest, "_cancel_run", mock_cancel_run)
#    monkeypatch.setattr(jaws_central.rest, "_cancel_transfer", mock_cancel_transfer)
#    with pytest.raises(Exception):
#        jaws_central.rest.cancel_run("user", 123)
#
#    """Check if no exception is raised in case run has a regular status"""
#    monkeypatch.setattr(jaws_central.rest, "_get_run", get_run_regular)
#    jaws_central.rest.cancel_run("user", 123)
#
#
def test_info(monkeypatch):
    mock_session = MockSession()
    mock_data = MockRunModel()
    run = runs.Run(mock_session, mock_data)

    expected_info_keys = [
        "id",
        "compute_site_id",
        "result",
        "status",
        "updated",
        "tag"
    ]
    expected_more_info_keys = [
        "submitted",
        "submission_id",
        "cromwell_run_id",
        "status_detail",
        "input_site_id",
        "upload_id",
        "download_id",
        "user_id",
        "wdl_file",
        "json_file",
    ]
    expected_info_keys_verbose = expected_info_keys + expected_more_info_keys

    complete_results = run.info(verbose=True)
    for key in complete_results:
        assert key in expected_info_keys_verbose
    for key in expected_info_keys_verbose:
        assert key in complete_results

    partial_results = run.info(verbose=False)
    for key in partial_results:
        assert key in expected_info_keys
    for key in expected_info_keys:
        assert key in partial_results


def test_inputs_manifest():
    mock_session = MockSession()
    mock_data = MockRunModel(manifest_json='{"foo": "bar"}')
    run = runs.Run(mock_session, mock_data)
    inputs_manifest = run.inputs_manifest()
    assert type(inputs_manifest) is dict
    assert inputs_manifest["foo"] == "bar"


def test_update_status(monkeypatch):
    def mock__update_status(self, status_from, status_to):
        assert status_from != status_to
        self.RUN_STATUS_WAS_UPDATED = True

    def mock__insert_run_log(self, status_from, status_to, timestamp, reason):
        assert status_from != status_to
        self.RUN_LOG_WAS_INSERTED = True

    monkeypatch.setattr(runs.Run, "_update_status", mock__update_status)
    monkeypatch.setattr(runs.Run, "_insert_run_log", mock__insert_run_log)

    mock_session = MockSession()
    mock_data = MockRunModel(status="queued")
    run = runs.Run(mock_session, mock_data)
    run.update_status("running")
    assert run.RUN_STATUS_WAS_UPDATED
    assert run.RUN_LOG_WAS_INSERTED
