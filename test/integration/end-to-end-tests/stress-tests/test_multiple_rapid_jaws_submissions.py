#!/usr/bin/env python

import json
import os
from datetime import datetime

import pytest
import submission_utils as util

# Note that in order for this test to do all the submissions at the same time, it MUST be run with
# pytest -n option with a value equal to the number of parameterized wdl/json combos.
# This means that increasing the number of wdl/json combinations will require updating -n value for
# the gitlab-ci.yml entries for this test


class TestMultipleRapidJawsSubmissions:
    FILE_DIR = '/global/cscratch1/sd/jaws/stress_tests'

    @pytest.mark.parametrize(
        "wdl, input_json, tag, num_of_submissions, check_tries, check_sleep",
        [
            (
                "../../../examples/jaws-alignment-example/alignment.wdl",
                "../../../examples/jaws-alignment-example/inputs.align.json",
                "alignment-stress-test",
                2, 360, 100,
            )
        ]
    )
    def test_rapid_submissions(
            self, env, site, wdl, input_json, tag, num_of_submissions, check_tries, check_sleep):

        # create a file named with the tag value plus a timestamp that will be used to record the
        # run_ids created by this test
        timestamp = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        filename = '%s/%s_%s.txt' % (TestMultipleRapidJawsSubmissions.FILE_DIR, tag, timestamp)

        # submit the wdl/json to JAWS num_of_submissions times and write the run ids to a file
        with open(filename, "w") as file:
            submitted = 0
            while submitted < num_of_submissions:
                jaws_output = util.submit_wdl(env, wdl, input_json, site, tag)
                run_id = str(jaws_output["run_id"])
                submitted += 1
                file.write('%s\n' % run_id)

        # change the file permissions to 775, so that it can be read by anyone later
        # not just the ci/cd test runner user, or the person who kicked tests off manually
        os.chmod(filename, 0o775)

        # read the file of run ids and check whether each run has completed
        with open(filename, "r") as file:
            lines = file.readlines()
            for line in lines:
                run_id = line.strip()
                util.wait_for_run_and_check_for_success(
                    run_id, env, check_tries, check_sleep
                )
