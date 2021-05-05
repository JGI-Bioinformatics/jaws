#!/usr/bin/env python
"""
These functions are to test the "testcases" from the "score_card" integration tests.
google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

This library of tests uses "fixtures" from conftest.py which should be located in the same directory. 
There is no need to import conftest.py as it is done automatically.

This test checks that expected resources are available. For example, are the large memory machines available (e.g. skylake) on cori?

Author: Jeff Froula <jlfroula@lbl.gov>
Updated: 02/23/21
"""

import pytest
import json
import submission_utils as util

#########################
###     Functions     ###
#########################
def test_skylake_250G(env, submit_skylake_250):
    """tests that the jaws history command has the run id in the stdout.
    This test will be skipped if the 'env' is not cori.
    """
    run_id = str(submit_skylake_250["run_id"])

    # Run: jaws status
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws status --verbose %s" % (env, run_id)
    (r, o, e) = util.run(cmd)

    status_info = json.loads(o)
    assert status_info["status"] == "download complete"
    assert status_info["result"] == "succeeded"


def test_skylake_500G(env, submit_skylake_500):
    """tests that the jaws history command has the run id in the stdout.
    This test will be skipped if the 'env' is not cori.
    """
    run_id = str(submit_skylake_500["run_id"])

    # Run: jaws status
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws status --verbose %s" % (env, run_id)
    (r, o, e) = util.run(cmd)

    status_info = json.loads(o)
    assert status_info["status"] == "download complete"
    assert status_info["result"] == "succeeded"
