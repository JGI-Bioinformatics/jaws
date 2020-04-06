import requests

from jaws_site import dispatch
from tests.conftest import CROMWELL_ID


def test_server_status_down(monkeypatch, server_status_down):
    monkeypatch.setattr(requests, "get", server_status_down)
    response = dispatch.server_status({})
    assert 'error' in response


def test_server_status_up(monkeypatch, server_status_up):
    monkeypatch.setattr(requests, "get", server_status_up)
    response = dispatch.server_status({})
    assert 'UP' in response['result']['Cromwell']


def test_get_metadata_from_a_run(monkeypatch, metadata_get):
    monkeypatch.setattr(requests, "get", metadata_get)
    params = {"cromwell_id": CROMWELL_ID}
    response = dispatch.run_metadata(params)
    assert 'result' in response


def test_get_task_status(monkeypatch, metadata_get):
    monkeypatch.setattr(requests, "get", metadata_get)
    params = {"cromwell_id": CROMWELL_ID}
    response = dispatch.task_status(params)
    print(response['result'])
    ordering = ['sc_test.do_prepare', 'sc_test.do_scatter',
                'sc_test.do_gather']
    result = response['result']

    # returns a tuple of job, status, start, end
    assert ordering[0] == result[0][0]
    assert ordering[1] == result[1][0]
    assert ordering[2] == result[2][0]


def test_get_task_ids(monkeypatch, metadata_get):
    monkeypatch.setattr(requests, "get", metadata_get)
    params = {"cromwell_id": CROMWELL_ID}
    response = dispatch.task_status(params)
    print(response['result'])

    results = response['result']

    all_tasks = ['sc_test.do_gather',
                 'sc_test.do_prepare',
                 'sc_test.do_scatter']

    assert len(results) == 3

    for task in results:
        job_name = task[0]
        assert job_name in all_tasks


def test_get_logs(monkeypatch, logs_get):
    monkeypatch.setattr(requests, "get", logs_get)
    params = {"cromwell_id": CROMWELL_ID}
    response = dispatch.run_logs(params)

    result = response['result']

    assert 'logs' in result


def test_cancel_run(monkeypatch, abort_post):
    monkeypatch.setattr(requests, "post", abort_post)
    params = {"cromwell_id": CROMWELL_ID}
    response = dispatch.cancel_run(params)

    assert "Cancelling" == response['result']['status']


def test_tail_files(log_file):
    logs, is_truncated = dispatch.tail(log_file)
    print(logs)
    assert is_truncated
    assert len(logs) == 1000
    assert logs[0] == 'this is line number 2000\n'
    assert logs[-1] == 'this is line number 2999\n'


def test_find_rc_failed_files(cromwell_run_dir):
    run_dir = cromwell_run_dir
    output = dispatch.find_failure_logs(run_dir)
    logs_with_failure = ["call-asm_1", "call-circularizeAssembly",
                         "call-filterHighGc"]
    print(output)
    for job in logs_with_failure:
        assert job in output
