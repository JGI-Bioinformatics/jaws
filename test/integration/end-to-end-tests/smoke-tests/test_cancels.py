#!/usr/bin/env python
"""
These functions are to test the "testcases" from the "score_card" integration tests.
google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

This library of tests uses "fixtures" from conftest.py which should be located in the same directory.
"""

import os
import pytest
import json
import time
import submission_utils as util


WDL = "/WDLs/fq_count.wdl"
INP = "/test-inputs/fq_count.json"
FINAL_STATE = "download complete"
VALID_STATES = [
    "queued",
    "running",
    "download complete"]
CHECK_TRIES = 20
CHECK_SLEEP = 10


@pytest.mark.parametrize("state", VALID_STATES)
def test_cancel(env, dir, site, state):
    wdl = dir + WDL
    input_json = dir + INP
    run_id = util.submit_wdl(env, wdl, input_json, site)["run_id"]

    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    status_cmd = source_cmd + "jaws status %s" % run_id
    cancel_cmd = source_cmd + "jaws cancel %s" % run_id

    ## get to the state we are testing
    floating_state = ""
    tries = 0
    while floating_state != state:
        r_sts, o_sts, e_sts = util.run(status_cmd)
        floating_state = json.loads(o_sts)["status"]
        if floating_state == FINAL_STATE:
            break

        if tries < CHECK_TRIES:
            time.sleep(CHECK_SLEEP)
        else:
            assert 0, "\n**Run took too long; last state: %s. (%s)" % (state, run_id)
            return
        tries += 1

    if floating_state != state and floating_state == FINAL_STATE:
        print("\n** This state was not observed: %s" % state)
        return

    ## submit cancel command to JAWS and gather output
    r_can, o_can, e_can = util.run(cancel_cmd)
    if state == FINAL_STATE:
        assert r_can == 1, "\n** Return code is not 1. (%s)" % (run_id)
        assert "error" in e_can, "\n** Error keyword absent. (%s)" % (run_id)
        return

    assert r_can == 0, "\n** Return Code is not 0. (%s)" % (run_id)
    assert "error" not in e_can, "\n** Error keyword present. (%s)" % (run_id)
    return