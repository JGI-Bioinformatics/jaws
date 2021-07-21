#!/usr/bin/env python

# These functions are to test the "score_card" unit tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
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

import sys
import os
import json
import re
import pytest
import time
import configparser
import glob
import submission_utils as util

check_tries = 50
check_sleep = 60

#########################
###     Functions     ###
#########################


def test_task_status(env, submit_subworkflow_alignment):
    """
    # task-status verifies all subworkflows task status was shown
    #
    #TASK_NAME  ATTEMPT CROMWELL_JOB_ID STATUS_FROM     STATUS_TO       TIMESTAMP       REASON  STATUS_DETAIL
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.shard      1 46806   running success 2021-02-08 20:53:55  The job completed successfully
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.bbmap_indexing       1       46807   running success 2021-02-08 20:53:58  The job completed successfully
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.alignment  1 46808   running success 2021-02-08 20:54:25  The job completed successfully
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.merge_bams 1 46809   running success 2021-02-08 20:54:37  The job completed successfully
    main_wdl.bam_stats                                1 46810   running success 2021-02-08 20:56:36  The job completed successfully
    """

    # time.sleep(120)  # wait some time before running task-status since there is some lag between
    #                 # when "jaws status" calls success and when "jaws task-status" calls success.
    run_id = submit_subworkflow_alignment["run_id"]
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws task-status %s | tail -n+2" % (
        env,
        run_id,
    )
    (r, o, e) = util.run(cmd)

    # put the table into a dictionary
    task_names = []
    status_to = []
    line_list = o.split("\n")
    line_list = list(filter(None, line_list))  # remove empty element
    for i in line_list:
        task_names.append(i.split("\t")[0])
        status_to.append(i.split("\t")[5])

    # print(f"TASK names: {status_to}")

    # check that the subworkflows tasks are in the list
    assert len(task_names) == 5

    # make sure all tasks completed with success
    assert len(list(filter(lambda x: (x == "success"), status_to))) == 5


def test_task_log(env, submit_subworkflow_alignment):
    """
    Test that all subworkflow tasks are represented by the task-log command

    #TASK_NAME  ATTEMPT CROMWELL_JOB_ID STATUS_FROM     STATUS_TO       TIMESTAMP       REASON
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.shard        1       46888   created ready   2021-02-09 22:04:14
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.shard        1       46888   ready   queued  2021-02-09 22:04:16
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.shard        1       46888   queued  pending 2021-02-09 22:04:17
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.shard        1       46888   pending running 2021-02-09 22:05:30
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.shard        1       46888   running success 2021-02-09 22:05:32
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.bbmap_indexing       1       46889   created ready   2021-02-09 22:04:20
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.bbmap_indexing       1       46889   ready   queued  2021-02-09 22:04:22
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.bbmap_indexing       1       46889   queued  pending 2021-02-09 22:04:23
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.bbmap_indexing       1       46889   pending running 2021-02-09 22:05:34
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.bbmap_indexing       1       46889   running success 2021-02-09 22:05:35
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.alignment    1       46890   created ready   2021-02-09 22:05:56
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.alignment    1       46890   ready   queued  2021-02-09 22:05:57
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.alignment    1       46890   queued  pending 2021-02-09 22:05:59
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.alignment    1       46890   pending running 2021-02-09 22:06:03
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.alignment    1       46890   running success 2021-02-09 22:06:06
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.merge_bams   1       46891   created ready   2021-02-09 22:06:11
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.merge_bams   1       46891   ready   queued  2021-02-09 22:06:12
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.merge_bams   1       46891   queued  pending 2021-02-09 22:06:14
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.merge_bams   1       46891   pending running 2021-02-09 22:06:18
    main_wdl.bbmap_shard_wf.bbmap_shard_wf.merge_bams   1       46891   running success 2021-02-09 22:06:19
    main_wdl.bam_stats  1       46892   created ready   2021-02-09 22:06:29
    main_wdl.bam_stats  1       46892   ready   queued  2021-02-09 22:06:31
    main_wdl.bam_stats  1       46892   queued  pending 2021-02-09 22:06:32
    main_wdl.bam_stats  1       46892   pending running 2021-02-09 22:07:30
    main_wdl.bam_stats  1       46892   running success 2021-02-09 22:07:31
    """

    # time.sleep(120)  # wait some time before running task-status since there is some lag between
    #                 # when "jaws status" calls success and when "jaws task-status" calls success.

    run_id = submit_subworkflow_alignment["run_id"]
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws task-log %s | tail -n+2" % (
        env,
        run_id,
    )
    (r, o, e) = util.run(cmd)

    # put the table into a dictionary
    task_names = []
    line_list = o.split("\n")
    line_list = list(filter(None, line_list))  # remove empty element
    for i in line_list:
        task_names.append(i.split("\t")[0])

    # check that the subworkflows tasks are in the list
    assert len(task_names) == 25


def test_for_raw_cromwell_files(env, submit_subworkflow_alignment):
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

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws status --verbose %s" % (env, run_id)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)
    outdir = data["output_dir"]

    cmd = (
        "find %s/call-bbmap_shard_wf/align.bbmap_shard_wf -name rc -exec cat {} \\; | grep -c 0"
        % (outdir)
    )
    (r, o, e) = util.run(cmd)

    # make sure all 4 of our "rc" files returned 0
    assert int(o.strip()) == 4


def test_saved_subwdl(env, submit_subworkflow_alignment):
    """
    subworkflow WDLs are saved in the user defined output dir

    """
    run_id = submit_subworkflow_alignment["run_id"]

    # need to get the submission_id from the status
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws status --verbose %s" % (env, run_id)
    (r, o, e) = util.run(cmd)
    assert not r
    data = json.loads(o)
    submission_id = data["submission_id"]
    outdir = data["output_dir"]

    zip_file = os.path.join(outdir, submission_id + ".zip")

    assert os.path.exists(zip_file)

    cmd = "unzip -l %s" % (zip_file)
    (r, o, e) = util.run(cmd)
    assert "alignment.wdl" in o


def test_subworkflow_metadata(env, submit_subworkflow_alignment):
    """
    metadata command also returns cromwell metadata for subworkflows
    """
    run_id = submit_subworkflow_alignment["run_id"]
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws metadata %s" % (env, run_id)
    (r, o, e) = util.run(cmd)
    meta_output = json.loads(o)

    # make sure metadata returns these calls
    expected = [
        "bbmap_shard_wf.alignment",
        "bbmap_shard_wf.bbmap_indexing",
        "bbmap_shard_wf.merge_bams",
        "bbmap_shard_wf.shard",
    ]
    calls = []
    for key in meta_output:
        for yek in meta_output[key]["calls"]:
            calls.append(yek)

    assert len([x for x in expected if x in calls]) == 4
