#!/usr/bin/env python

import json
import pytest
import submission_utils as util


class TestRunSuccess:
    check_tries = 360  # try this many times when waiting for a JAWS run to complete.
    check_sleep = 100  # wait for this amount of time between tries.

    @staticmethod
    def run_success(env, site, wdl, input_json):

        jaws_output = util.submit_wdl(env, wdl, input_json, site)
        run_id = str(jaws_output["run_id"])

        util.wait_for_run(
            run_id, env, TestRunSuccess.check_tries, TestRunSuccess.check_sleep
        )

        cmd = "source ~/jaws-%s.sh > /dev/null && jaws status %s" % (env, run_id)

        (rc, stdout, stderr) = util.run(cmd)
        # print("status cmd:", cmd)
        # print("rc: ", rc)
        # print("stderr: ", stderr)
        # print("stdout: ", stdout)

        status_info = json.loads(stdout)
        assert status_info["status"] == "download complete",  \
            "\n**Run %s took too long - last state seen: %s" % (run_id, status_info["status"])
        assert status_info["result"] == "succeeded", \
            "\n**Run %s did not succeed - status was: %s" % (run_id, status_info["result"])


class TestMultipleRapidJawsSubmissions(TestRunSuccess):
    @pytest.mark.parametrize(
        "wdl, input_json",
        (
            (
                "../../../examples/jaws-alignment-example/alignment.wdl",
                "../../../examples/jaws-alignment-example/inputs.align.json",
            ),
        ),
    )
    def test_tutorial_success(self, clone_tutorials_repo, env, site, wdl, input_json):
        # run the test against all the wdl/json in the @pytest.mark.parametrize
        self.run_success(env, site, wdl, input_json)
