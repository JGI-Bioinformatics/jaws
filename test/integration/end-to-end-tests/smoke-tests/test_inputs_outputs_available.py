#!/usr/bin/env python

# These functions are to test the "testcases" from the "score_card" integration tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n  # noqa
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.
# There is no need to import conftest.py as it is done automatically.

import os
import uuid
import shutil
import submission_utils as util

# set variables specific for this series of tests
check_tries = 100
check_sleep = 30

#####################
#     Functions     #
#####################


def test_jaws_get(submit_fq_count_wdl):
    """
    1) I will check that the input WDL and json file are saved in the output dir.
       (70f82d8f-352f-48ce-a21d-3e4ede4daef3.orig.json  70f82d8f-352f-48ce-a21d-3e4ede4daef3.wdl)

    2) I will also check that the raw cromwell file structure was created

    3) I will check that the expected files were created in the execution dir (i.e. stdout and stderr)
        This cromwell file structure should exist
        fq_count_out/call-count_seqs/execution/
            num_seqs.txt  rc  script  script.submit  stderr  stderr.submit      stdout  stdout.submit
    """
    input_wdl = "main.wdl"
    input_json = "inputs.json"
    outdir = str(uuid.uuid4())

    run_id = str(submit_fq_count_wdl["run_id"])
    util.wait_for_run(run_id, check_tries, check_sleep)

    cmd = "jaws get --quiet --complete %s %s" % (run_id, outdir)
    (r, o, e) = util.run(cmd)
    assert r == 0

    # check that we have the initial WDL saved to the outdir
    # using the full path (outdir and input_wdl), we are essentially testing that the outdir
    # was correct and that the wdl file got created.

    # verify it is a valid wdl
    with open(os.path.join(outdir, input_wdl)) as fh:
        # check of the workflow title is in the top 3 lines
        head = [next(fh) for x in range(3)]
        head = ' '.join(head)
        if "workflow fq_count" not in head:
            assert 0, "This does not look like a valid workflow"

    # check that we have a valid inputs json
    with open(os.path.join(outdir, input_json)) as fh:
        expected = '"fq_count.fastq_file":'
        if expected not in fh.read():
            assert 0, "This does not look like a valid inputs json file"

    expected_files = [
        "num_seqs.txt",
        "rc",
        "script",
        "script.submit",
        "stderr",
        "stderr.submit",
        "stdout",
        "stdout.submit",
    ]
    for file in expected_files:
        if not os.path.exists(os.path.join(outdir, "call-count_seqs/execution/", file)):
            assert (
                0
            ), f"expected result file not found in: {os.path.join(outdir,'call-count_seqs/execution/',file)}"

    try:
        shutil.rmtree(outdir)
    except OSError as error:
        print(f"Error: {outdir}: {error}")
