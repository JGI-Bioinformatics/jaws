#!/usr/bin/env python

# These functions are to test the "score_card" unit tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n # noqa
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451
#
# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.
#
# The following tests will make sure JAWS handles subworkflows correctly.
# 1) task-status verifies all subworkflows task status was shown
# 2) task-log verifies all subworkflows task status was shown
# 3) raw cromwell subworkflow files are returned to user defined output dir
# 4) subworkflow WDLs are saved in the user defined output dir
# 5) metadata command also returns cromwell metadata for subworkflows

import os
import json
import shutil
import uuid
import submission_utils as util

check_tries = 50
check_sleep = 60

#####################
#     Functions     #
#####################


def test_for_raw_cromwell_files(submit_subworkflow_alignment):
    """
    test that raw cromwell subworkflow files are returned to user defined output dir.

    These files should exist OUTDIR/call-bbmap_shard_wf/<cromwell-hash>/
        shard
        bbmap_indexing
        alignment
        merge_bams

        out/call-bbmap_shard_wf/align.bbmap_shard_wf/c9f67f71-1acd-4d8c-8879-14b7b1a28b54/call-shard/execution/rc
    """
    run_id = submit_subworkflow_alignment["run_id"]

    outdir = str(uuid.uuid4())

    cmd = "jaws get --quiet --complete %s %s" % (run_id, outdir)
    (r, o, e) = util.run(cmd)
    assert not r

    cmd = (
        "find %s/call-bbmap_shard_wf/align.bbmap_shard_wf -name rc -exec cat {} \\; | grep -c 0"
        % (outdir)
    )
    (r, o, e) = util.run(cmd)

    # make sure all 4 of our "rc" files returned 0
    assert int(o.strip()) == 4

    try:
        shutil.rmtree(outdir)
    except OSError as error:
        print(f"Error: {outdir}: {error}")

@pytest.mark.xfail(reason="We have a ticket for this")
def test_saved_subwdl(submit_subworkflow_alignment):
    """
    subworkflow WDLs are saved in the user defined output dir

    """
    run_id = submit_subworkflow_alignment["run_id"]
    outdir = str(uuid.uuid4())
    cmd = "jaws get --quiet --complete %s %s" % (run_id, outdir)
    (r, o, e) = util.run(cmd)
    assert not r

    zip_file = os.path.join(outdir, f"subworkflows.zip")
    assert os.path.exists(zip_file)

    cmd = "unzip -l %s" % (zip_file)
    (r, o, e) = util.run(cmd)
    assert "alignment.wdl" in o

    try:
        shutil.rmtree(outdir)
    except OSError as error:
        print(f"Error: {outdir}: {error}")


def test_subworkflow_metadata(submit_subworkflow_alignment):
    """
    metadata command also returns cromwell metadata for subworkflows
    """
    run_id = submit_subworkflow_alignment["run_id"]
    cmd = "jaws metadata %s" % (run_id)
    (r, o, e) = util.run(cmd)
    meta_output = json.loads(o)

    # make sure metadata returns these calls
    expected = [
        "bbmap_shard_wf.alignment",
        "bbmap_shard_wf.bbmap_indexing",
        "bbmap_shard_wf.merge_bams",
        "bbmap_shard_wf.shard",
    ]
    calls = meta_output["calls"]["main_wdl.bbmap_shard_wf"][0]["subWorkflowMetadata"][
        "calls"
    ].keys()
    assert len([x for x in expected if x in calls]) == 4
