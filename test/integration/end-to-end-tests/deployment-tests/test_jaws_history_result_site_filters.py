#!/usr/bin/env python
import pytest
import json
import submission_utils as util


def test_jaws_history_site_filter(env, site):
    """
    jaws history --site [CORI, JGI, CASCADE]
    """
    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws history --site %s" % (env, site)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
      for d in data:
          assert d["site_id"] == site
      else:
          print("No output from jaws history command")
          assert 1 == 1


def test_jaws_history_result_filter(env):
    """
    jaws history --result [succeeded, failed]
    Checking the output only with "succeeded" and "failed"
    """
    for k in ["succeeded", "failed"]:
        cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws history --result %s" % (env, k)
        (r, o, e) = util.run(cmd)
        data = json.loads(o)

        if data:
            for d in data:
              assert d["result"] == k
        else:
            print("No output from jaws history command")
            assert 1 == 1

