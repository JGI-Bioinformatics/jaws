import os
import json
import pytest
import time
from subprocess import Popen, PIPE

def run(cmd):
    output = Popen(cmd, stdout=PIPE,
             stderr=PIPE, shell=True,
             universal_newlines=True)

    stdout,stderr=output.communicate()
    rc=output.returncode

    return rc,stdout,stderr


def submit_wdl(env, wdl, input_json, outdir, site):
    """
    This is a fixture that will submit a wdl for all functions to use.
    This function returns the output of a wdl submission.
    """
    if os.path.exists(outdir):
        cmd = "rm -rf %s" % outdir
        print(cmd)
        (rc, stdout, stderr) = run(cmd)
        if rc > 0:
            os.exit("Failed to remove old output directory %s" % outdir)

    cmd = "source ~/jaws-%s.sh && jaws run submit %s %s %s %s" % (env, wdl, input_json, outdir, site)
    print(cmd)
    (rc, stdout, stderr) = run(cmd)
    if rc > 0:
        pytest.exit("stderr: %s    stdout: " % stderr, stdout)

    assert rc == 0
    data = json.loads(stdout)

    # uncomment for testing
    # data={
    #    "output_dir": "fq_count_out",
    #    "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
    #    "run_id": 16405,
    #    "site_id": "CORI",
    #    "status": "uploading",
    #    "submission_id": "47a555c7-07a6-442c-a2f1-d0319f2e3008",
    #    "upload_task_id": "444ac0b8-60f0-11eb-9905-0aa9ddbe2755"
    # }
    return data