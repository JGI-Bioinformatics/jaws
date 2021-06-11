#!/usr/bin/env python

# These functions are to test the "testcases" from the "score_card" integration tests.
# google doc: https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n
# google sheet: https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451

# This library of tests uses "fixtures" from conftest.py which should be located in the same directory.

import os
import pytest
import json
import time
import submission_utils as util
import tempfile
import logging
import socket
import re

#log = logging.getLogger('test_cache.log',level=INFO)

check_tries = 100
check_sleep = 30


#########################
###     Functions     ###
#########################

def get_submission_id(env,cmd):
    (r, o, e) = util.run(cmd)
    if r != 0: 
        print("jaws submit command failed %s" % cmd)
        assert 0

#    o="""
#    test_cache.py::test_cache MY OUT: {
#    "max_ram_gb": 10,
#    "run_id": 7615,
#    "site_id": "CORI",
#    "status": "uploading",
#    "tag": null
#}"""

    # The stdout from a submission looks like
    #
    # test_cache.py::test_cache: {
    # "max_ram_gb": 10,
    # "run_id": 7615,
    # "site_id": "CORI",
    # "status": "uploading",
    # "tag": null
    # }
    x = re.search(r"{.*}",o,re.DOTALL)
    run_id = 0
    if x:
        json_only = x.group(0)
        data = json.loads(json_only)
        run_id = data["run_id"]
    else:
        pytest.exit(f"there was no json in the output for the submit command: {cmd}")

    util.wait_for_run(run_id, env, check_tries, check_sleep)

    cmd = (
        "source ~/jaws-%s.sh > /dev/null && jaws metadata %s" 
        % (env, run_id)
    )
    (r, o, e) = util.run(cmd)
    if r != 0: 
        print("jaws metadata command failed %s" % cmd)
        assert 0
    data = json.loads(o)
    return data,run_id

def test_cache(env,site,dir):
    """make sure cromwell's cache'ing feature works and can be turned off with jaws submit --no-cache."""

    wdl = dir + '/WDLs/fq_count.wdl'
    input_json_file = dir + '/test_cache.json'

    # By creating a unique name for the inputs json file, we should force Cromwell to not use a cached result.
    # We should see the callCaching: hit equal to false
    # i.e. data['submission_id']['calls']['fq_count.count_seqs'][0]['callCaching']['hit'] == false
    randstring = tempfile.NamedTemporaryFile(prefix='sample.',dir='.')
    temp_json_path = "%s.fastq" % randstring.name
    myhost = socket.gethostname()
    if 'cori' in myhost:
        os.symlink('/global/cfs/projectdirs/jaws/test/tutorial_test_data/sample.fastq', temp_json_path)
    elif 'lrc' in myhost:
        os.symlink('/global/home/groups-sw/lr_jgicloud/examples_data/sample.fastq', temp_json_path)
    else:
        sys.exit('Error: hostname unidentified')

    content="{\"fq_count.fastq_file\": \"%s\"}" % temp_json_path
    with open(input_json_file,'w') as mj:
        mj.write(content)
        
    cmd = (
        "source ~/jaws-%s.sh > /dev/null && jaws submit %s %s %s" 
        % (env, wdl, input_json_file, site)
    )

    (data,run_id) = get_submission_id(env,cmd)
    submission_id = list(data.keys())[0]

    # we should not see caching with this submission
    assert data[submission_id]['calls']['fq_count.count_seqs'][0]['callCaching']['hit'] is False, "This run \"%s\" should not be cached but was." % run_id


    ### ------------------------------------------------- ###    
    # Running again with the same input should use the cached result.
    cmd = (
        "source ~/jaws-%s.sh > /dev/null && jaws submit %s %s %s" 
        % (env, wdl, input_json_file, site)
    )

    (data,run_id) = get_submission_id(env,cmd)
    submission_id = list(data.keys())[0]

    # we should not see caching with this submission
    assert data[submission_id]['calls']['fq_count.count_seqs'][0]['callCaching']['hit'] is True, "This run \"%s\" should be cached but was not." % run_id


    ### ------------------------------------------------- ###    
    # Test that --no-cache works with same input
    cmd = (
        "source ~/jaws-%s.sh > /dev/null && jaws submit --no-cache %s %s %s" 
        % (env, wdl, input_json_file, site)
    )

    (data,run_id) = get_submission_id(env,cmd)
    submission_id = list(data.keys())[0]

    # we should not see caching with this submission
    assert data[submission_id]['calls']['fq_count.count_seqs'][0]['callCaching']['allowResultReuse'] is False, "This run \"%s\" should not be cached (--no-cache used) but was." % run_id
