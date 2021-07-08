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

# Also note that this test takes a --site-list parameter, rather than --site (so that submissions
# to all sites takes place simultaneously).  Here is an example of how to run this test from
# the command line:
# > pytest -n 4 --capture=no --verbose --dir . --site-list "cori, jgi" --env staging
#   stress-tests/test_multiple_rapid_jaws_submissions.py

class TestMultipleRapidJawsSubmissions:
    FILE_DIR = '/global/cscratch1/sd/jaws/stress_tests'

    # note that the num_of_submissions is the number of submissions to EACH site in the
    # site-list param, so, if there are 2 sites in the list the total number of submissions to
    # jaws will be 2 * (the sum of all the num_of_submissions in the list)
    @pytest.mark.parametrize(
        "wdl, input_json, tag, num_of_submissions, check_tries, check_sleep",
        [
            (
                "../../../examples/bfoster_meta_assem/jgi_meta.jaws.wdl",
                "../../../examples/bfoster_meta_assem/inputs.json",
                "bfoster-small-job-stress-test",
                30, 360, 120,
            ),
            (
                "../../../examples/leo_dapseq/Azospirillum_brasilense.wdl",
                "../../../examples/leo_dapseq/shortened.json",
                "leo-small-job-stress-test",
                10, 360, 120,
            ),
            (
                "../../../examples/bfoster_meta_assem/jgi_meta.jaws.wdl",
                "../../../examples/bfoster_meta_assem/big_inputs.json",
                "bfoster-big-job-stress-test",
                5, 360, 120,
            ),
            (
                "../../../examples/leo_dapseq/Azospirillum_brasilense.wdl",
                "../../../examples/leo_dapseq/shortened-100.json",
                "leo-medium-job-stress-test",
                5, 360, 120,
            )
        ]
    )
    def test_rapid_submissions(
            self, env, site_list, wdl, input_json, tag, num_of_submissions, check_tries, check_sleep):

        # site_list will be comma seperated string, like "cori, jgi"
        # make it a python list with no blanks
        sites = [x.strip() for x in site_list.split(',')]

        # create a file named with the tag value plus a timestamp that will be used to record the
        # run_ids created by this test
        timestamp = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        filename = '%s/%s_%s.txt' % (TestMultipleRapidJawsSubmissions.FILE_DIR, tag, timestamp)

        # submit the wdl/json to JAWS num_of_submissions times and write the run ids to a file
        with open(filename, "w") as file:
            # submit to each site in site_list
            submitted = 0
            while submitted < num_of_submissions:
                for site in sites:
                    jaws_output = util.submit_wdl(env, wdl, input_json, site, tag)
                    run_id = str(jaws_output["run_id"])
                    file.write('%s\n' % run_id)

                submitted += 1

        # change the file permissions to 775, so that it can be read by anyone later
        # not just the ci/cd test runner user or the person who kicked tests off manually
        os.chmod(filename, 0o775)

        # read the file of run ids and check whether each run has completed
        with open(filename, "r") as file:
            lines = file.readlines()
            for line in lines:
                run_id = line.strip()
                util.wait_for_run_and_check_for_success(
                    run_id, env, check_tries, check_sleep
                )
