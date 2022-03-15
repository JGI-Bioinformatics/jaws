from dataclasses import dataclass
from datetime import datetime
from jaws_site import tasks, runs, runs_es, cromwell

this_date = datetime.today()


def mock_task_summary(*args, **kwargs):
    return [
        [
            'task_abcd',      # task_name
            'cromwell_abcd',  # cromwell_job_id
            False,            # cached
            'success',        # result
            this_date,        # queued
            '01:00:00',       # queue-wait
            '02:00:00',       # run-time
            '03:00:00',       # max-time
        ],
        [
            'task_efgh',      # task_name
            'cromwell_efgh',  # cromwell_job_id
            True,             # cached
            'success',        # result
            this_date,        # queued
            '04:00:00',       # queue-wait
            '05:00:00',       # run-time
            '06:00:00',       # max-time
        ],
    ]


def mock_task_status(*args, **kwargs):
    return [
        [
            'task_abcd',      # task_name
            'cromwell_abcd',  # cromwell_job_id
            'cached_abcd',    # cached
            'success',        # status
            this_date,        # timestamp
            'reason_abcd',    # reason
        ],
        [
            'task_efgh',      # task_name
            'cromwell_efgh',  # cromwell_job_id
            'cached_efgh',    # cached
            'success',        # status
            this_date,        # timestamp
            'reason_efgh',    # reason
        ],
    ]


def mock_metadata(*args, **kwargs):
    return {
        'workflow_name': 'test'
    }


def mock_cromwell_metadata(*args, **kwargs):
    @dataclass
    class Metadata:
        data = {'workflowName': 'test'}
    return Metadata()


def test_task_summary(mock_db_session, monkeypatch):
    monkeypatch.setattr(tasks.TaskLog, 'task_summary', mock_task_summary)

    exp_results = {
        'task_abcd':
            {
                'cached': False,
                'cromwell_id': 'cromwell_abcd',
                'maxtime': '03:00:00',
                'queue_wait': '01:00:00',
                'queued': this_date,
                'result': 'success',
                'runtime': '02:00:00'},
        'task_efgh':
            {
                'cached': True,
                'cromwell_id': 'cromwell_efgh',
                'maxtime': '06:00:00',
                'queue_wait': '04:00:00',
                'queued': this_date,
                'result': 'success',
                'runtime': '05:00:00'
            }
        }

    obs_results = runs_es.RunES(mock_db_session, run_id=123).task_summary()
    assert obs_results == exp_results


def test_task_status(mock_db_session, monkeypatch):
    monkeypatch.setattr(tasks.TaskLog, 'task_status', mock_task_status)

    exp_results = {
        'task_abcd': {
            'cromwell_id': 'cromwell_abcd',
            'reason': 'reason_abcd',
            'status': 'success',
            'timestamp': this_date,
        },
        'task_efgh': {
            'cromwell_id': 'cromwell_efgh',
            'reason': 'reason_efgh',
            'status': 'success',
            'timestamp': this_date,
        }
    }

    obs_results = runs_es.RunES(mock_db_session, run_id=123).task_status()
    assert obs_results == exp_results


def test_cromwell_metadata(mock_db_session, monkeypatch):
    exp_results = {
        'workflowName': 'test'
    }

    def mock_metadata(*args, **kwargs):
        return exp_results

    mock_db_session.output([{'cromwell_run_id': 'abc'}])
    monkeypatch.setattr(cromwell.Cromwell, 'get_metadata', mock_cromwell_metadata)

    obs_results = runs_es.RunES(mock_db_session, run_id=123).cromwell_metadata()
    assert obs_results == exp_results


def test_create_doc(mock_db_session, monkeypatch):
    monkeypatch.setattr(tasks.TaskLog, 'task_summary', mock_task_summary)
    monkeypatch.setattr(tasks.TaskLog, 'task_status', mock_task_status)
    monkeypatch.setattr(runs.Run, 'metadata', mock_metadata)
    monkeypatch.setattr(runs_es.RunES, 'cromwell_metadata', mock_cromwell_metadata)

    session_result = [
        {
            'id': 123,
            'run_id': 123,
            'user_id': 'John Doe',
            'email': 'johndoe@lbl.gov',
            'submitted': this_date,
            'updated': this_date,
            "status": "download complete",
            "result": "succeeded",
            "cromwell_run_id": "abcd"
        }
    ]
    mock_db_session.output(session_result, repeat=True)

    exp_results = {
        'email': 'johndoe@lbl.gov',
        'result': None,
        'run_id': 123,
        'site_id': 'EAGLE',
        'status': 'download complete',
        'status_detail': '',
        'submitted': this_date.strftime("%Y-%m-%d %H:%M:%S"),
        'tasks': [
            {
                'cached': False,
                'cromwell_id': 'cromwell_abcd',
                'maxtime': '03:00:00',
                'name': 'task_abcd',
                'queue_wait': '01:00:00',
                'queued': this_date,
                'reason': 'reason_abcd',
                'result': 'success',
                'runtime': '02:00:00',
                'status': 'success',
                'timestamp': this_date,
            },
            {
                'cached': True,
                'cromwell_id': 'cromwell_efgh',
                'maxtime': '06:00:00',
                'name': 'task_efgh',
                'queue_wait': '04:00:00',
                'queued': this_date,
                'reason': 'reason_efgh',
                'result': 'success',
                'runtime': '05:00:00',
                'status': 'success',
                'timestamp': this_date,
            }
        ],
        'updated': this_date.strftime("%Y-%m-%d %H:%M:%S"),
        'user_id': 'John Doe',
        'workflow_name': 'test'
    }
    obs_results = runs_es.RunES(mock_db_session, run_id=123).create_doc()
    if obs_results['site_id']:
        obs_results['site_id'] = obs_results['site_id'].upper()
    assert obs_results == exp_results


def test_send_rpc_run_metadata(mock_rpc_request):
    payload = {
        'result': 'okay'
    }
    exp_results = (None, 0)
    obs_results = runs_es.send_rpc_run_metadata(mock_rpc_request, payload)
    assert obs_results == exp_results
