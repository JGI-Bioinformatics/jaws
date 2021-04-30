#!/usr/bin/env python
import pytest
import json
import submission_utils as util


def test_jaws_queue_site_filter(env, site):
    """
    jaws queue --site [CORI, JGI, CASCADE]
    """
    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws queue --site %s" % (env, site)
    (r, o, e) = util.run(cmd)
    data = json.loads(o)

    if data:
        for d in data:
            assert d["site_id"] == site
    else:
        print("No output from jaws queue command")
        assert 1 == 1

