#!/usr/bin/env python

import json
import pytest
import submission_utils as util

# Note that in order for this test to do all the submissions at the same time, it MUST be run with
# pytest -n option with a value equal to the number of parameterized wdl/json combos.
# This means that increasing the number of wdl/json combinations will require updating -n value for
# the gitlab-ci.yml entries for this test


class TestMultipleRapidJawsSubmissions:

    @pytest.mark.parametrize(
        "wdl, input_json", "tag", "check_tries", "check_sleep",
        [
            (
                "../../../examples/jaws-alignment-example/alignment.wdl",
                "../../../examples/jaws-alignment-example/inputs.align.json",
                "alignment-stress-test",
                360, 100,
            )
        ]
    )
    def test_rapid_submissions(self, env, site, wdl, input_json, tag, check_tries, check_sleep):
        # run the test against all the wdl/json in the @pytest.mark.parametrize

        jaws_output = util.submit_wdl(env, wdl, input_json, site, tag)
        run_id = str(jaws_output["run_id"])

        util.wait_for_run(
            run_id, env, check_tries, check_sleep
        )

        cmd = "source ~/jaws-%s.sh > /dev/null && jaws status %s" % (env, run_id)

        (rc, stdout, stderr) = util.run(cmd)

        status_info = json.loads(stdout)
        assert status_info["status"] == "download complete", \
            "\n**Run %s took too long - last state seen: %s" % (run_id, status_info["status"])
        assert status_info["result"] == "succeeded", \
            "\n**Run %s did not succeed - status was: %s" % (run_id, status_info["result"])
