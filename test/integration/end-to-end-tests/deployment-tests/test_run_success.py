#!/usr/bin/env python

import json
import pytest
import submission_utils as util

check_tries = 100  # try this many times when waiting for a JAWS run to complete.
check_sleep = 30  # wait for this amount of time between tries.

@pytest.mark.parametrize(
    'wdl, input_json',
    (
            ('../../../../examples/bfoster_meta_assem/jgi_meta.jaws.wdl',
             '../../../../examples/bfoster_meta_assem/inputs.json'),
            ('./WDLs/jaws-alignment-example/main.wdl',
             './WDLs/jaws-alignment-example/inputs.json' )
    )
)
def test_run_success(env, site, wdl, input_json):

    jaws_output = util.submit_wdl(env, wdl, input_json, site)
    run_id = str(jaws_output['run_id'])

    util.wait_for_run(env, run_id, check_tries, check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run status %s" % (env, run_id)
    (rc, stdout, stderr) = util.run(cmd)
    print("%s\n%s\n%s\n%s", cmd, rc, stdout, stderr)

    status_info = json.loads(stdout)
    assert status_info['status'] == "download complete"
    assert status_info['result'] == "succeeded"