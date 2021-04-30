#!/usr/bin/env python
import pytest
import json
import submission_utils as util

check_tries=20
check_sleep=300

def test_jaws_history_site_filter(env, site, submit_fq_count_wdl):
    """
    jaws history --site [CORI, JGI, CASCADE]
    """
    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws history --site %s" % (env, site)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["site_id"].lower() == site
    else:
        pytest.exit(f"no runs were found in the history for site: {site}")


def test_jaws_history_result_filter_succeeded(env,submit_fq_count_wdl):
    """
    jaws history --result [succeeded, failed]
    Checking the output only with "succeeded" and "failed"
    """

    run_id = submit_fq_count_wdl['run_id']
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws history --result succeeded" % (env)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["result"] == 'succeeded'
    else:
        pytest.exit(f"no runs were found in the history that have result: succeeded")


def test_jaws_history_result_filter_failed(env,submit_bad_task):
    """
    jaws history --result failed
    Checking the output only with "succeeded" and "failed"
    """
    run_id = submit_bad_task['run_id']
    util.wait_for_run(run_id,env,check_tries,check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws history --result failed" % (env)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["result"] == 'failed'
    else:
        pytest.exit(f"no runs were found in the history that have result: failed")

