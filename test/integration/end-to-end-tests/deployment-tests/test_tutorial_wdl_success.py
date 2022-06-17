#!/usr/bin/env python
# The purpose of this test is more to make sure the tutorials are working than to
# test that JAWS is working.  However, there are some tests in here that are usefull
# to test JAWS components like making sure we can access the /refdata (referencing_db_and_shifter)

import pytest
import submission_utils as util


@pytest.mark.parametrize(
    "wdl, input_json",
    (
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
            "./jaws-tutorial-examples/copy-refdata-as-inputs/refdata.wdl",
            "./jaws-tutorial-examples/copy-refdata-as-inputs/refdata.json",
        ),
        (
            "./jaws-tutorial-examples/referencing_db_and_shifter/test.wdl",
            "./jaws-tutorial-examples/referencing_db_and_shifter/inputs.json",
        ),
        (
            "./jaws-tutorial-examples/scatter_gather_example/array_scatter.wdl",
            "./jaws-tutorial-examples/scatter_gather_example/array_scatter.json",
        ),
        (
            "./jaws-tutorial-examples/scatter_gather_example/map_scatter.wdl",
            "./jaws-tutorial-examples/scatter_gather_example/map_scatter.json",
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
def test_tutorial_success(clone_tutorials_repo, site, wdl, input_json):
    # run the test against all the wdl/json in the @pytest.mark.parametrize
    util.run_success(site, wdl, input_json)
