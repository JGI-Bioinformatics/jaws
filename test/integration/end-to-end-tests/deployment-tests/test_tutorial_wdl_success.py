#!/usr/bin/env python

import json
import pytest
import submission_utils as util


class TestRunSuccess:
    check_tries = 360  # try this many times when waiting for a JAWS run to complete.
    check_sleep = 60  # wait for this amount of time between tries.

    @staticmethod
    def run_success(env, site, wdl, input_json):

        jaws_output = util.submit_wdl(env, wdl, input_json, site)
        run_id = str(jaws_output["run_id"])

        util.wait_for_run_and_check_for_success(
            run_id, env, TestRunSuccess.check_tries, TestRunSuccess.check_sleep
        )


class TestTutorialSuccess(TestRunSuccess):
    @pytest.mark.parametrize(
        "wdl, input_json",
        (
            (
                "./jaws-tutorial-examples/fq_count/fq_count.wdl",
                "./jaws-tutorial-examples/fq_count/fq_count.json",
            ),
            (
                "./jaws-tutorial-examples/jaws-alignment-example/main.wdl",
                "./jaws-tutorial-examples/jaws-alignment-example/inputs.json",
            ),
            (
                "./jaws-tutorial-examples/jaws-sharding/shard_test.wdl",
                "./jaws-tutorial-examples/jaws-sharding/shard.json",
            ),
            (
                "./jaws-tutorial-examples/quickstart/align.wdl",
                "./jaws-tutorial-examples/quickstart/inputs.json",
            ),
            (
                "./jaws-tutorial-examples/referencing_db_and_shifter/test.wdl",
                "./jaws-tutorial-examples/referencing_db_and_shifter/inputs.json",
            ),
            (
                "./jaws-tutorial-examples/scatter_gather_example/test.wdl",
                "./jaws-tutorial-examples/scatter_gather_example/input.json",
            ),
            (
                "./jaws-tutorial-examples/subworkflow/main.wdl",
                "./jaws-tutorial-examples/subworkflow/inputs.json",
            ),
            (
                "./jaws-tutorial-examples/subworkflows_and_conditionals/main.wdl",
                "./jaws-tutorial-examples/subworkflows_and_conditionals/inputs.json",
            ),
        ),
    )
    def test_tutorial_success(self, clone_tutorials_repo, env, site, wdl, input_json):
        # run the test against all the wdl/json in the @pytest.mark.parametrize
        self.run_success(env, site, wdl, input_json)
