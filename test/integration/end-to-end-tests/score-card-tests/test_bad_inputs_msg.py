#!/usr/bin/env python
"""This will submit a WDL to JAWS and wait for completion, 
   run tests and finally create a yaml file for the autoqc GUI."""

import sys
import os
import argparse
import json
import logging
import parsing_functions as pf

check_tries = 50  # try this many times when waiting for a JAWS run to complete.
check_sleep = 30  # wait for this amount of time between tries.

def test_json_file_does_not_exist(env):
    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    wdl = "~/jaws/examples/jaws-alignment-example/main.wdl"
    out_dir = pf.timestamp_dir("./out/bad_inputs")

    # Submit job with path to json file that does not exist
    # Can't use pf.submit_one_run_to_env here because it exits if submission not successful
    json = "./FileDoesNotExist.json"

    submit_cmd = "jaws run submit %s %s %s %s" % (wdl, json, out_dir, args.site)
    cmd = source_cmd + submit_cmd
    (o, e, r) = pf.submit_cmd(cmd)
    # print("cmd: %s\nout: %s\nerror: %s", cmd, o, e)

    # check for the correct error message
    assert "No such file or directory:" in e


def test_json_bad_path_to_input_file_msg(env):
    source_cmd = "source ~/jaws-%s.sh > /dev/null && " % env
    wdl = "~/jaws/examples/jaws-alignment-example/main.wdl"
    out_dir = pf.timestamp_dir("./out/bad_inputs")

    # Submit job with json contains path to a non-existent input file
    # Can't use pf.submit_one_run_to_env here because it exits if submission not successful
    json = "../test-inputs/bad_path_inputs.json"
    submit_cmd = "jaws run submit %s %s %s %s" % (wdl, json, out_dir, args.site)
    cmd = source_cmd + submit_cmd
    (o, e, r) = pf.submit_cmd(cmd)
    # print("cmd: %s\nout: %s\nerror: %s", cmd, o, e)

    # check for the correct error message
    assert "File(s) not accessible:" in e

    # TODO add more tests to cover the other scorecard scenarios that I said this test would cover
