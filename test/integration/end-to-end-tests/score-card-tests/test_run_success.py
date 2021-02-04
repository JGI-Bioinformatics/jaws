#!/usr/bin/env python
"""This will submit a WDL to JAWS and wait for completions.
   It will then submit a test cmd and create an analysis.yaml file."""

import sys
import os
import argparse
import json
import logging
import parsing_functions as pf

check_tries = 100  # try this many times when waiting for a JAWS run to complete.
check_sleep = 30  # wait for this amount of time between tries.

def test_run_success(env, submit_wdl_and_wait):

    run_id = submit_wdl_and_wait['run_id']

    # Run: jaws run status
    cmd = "source ~/jaws-%s.sh > /dev/null && jaws run status %s" % (env, run_id)
    print(f"{cmd}\n")
    (o, e, r) = pf.submit_cmd(cmd)
    print("%s\n%s\n%s", cmd, o, e)

    status_info = json.loads(o)
    assert status_info['status'] == "download complete"
    assert status_info['result'] == "succeeded"