import os,sys
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
    if rc:
        sys.stderr.write("Error: This command returned an error code greater than 0: \n%s\n%s" % (cmd,stderr))

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

    # the pipe > /dev/null 2>&1 is needed below because otherwise the info printed from the
    # activation command causes an error when we try to do json load later

    cmd = "source ~/jaws-%s.sh > /dev/null 2>&1 && jaws run submit %s %s %s %s" % (env, wdl, input_json, outdir, site)
    (rc, stdout, stderr) = run(cmd)
    if rc > 0:
        if stderr:
            pytest.exit("stderr: %s" % stderr)
        if stdout:
            pytest.exit("stdout: %s" % stdout)
        else:
            pytest.exit("The return code from the command below returned something > 0. There was no stderr accompaning this failure. \n%s" % cmd)

    assert rc == 0
    data = json.loads(stdout)


    """
    # uncomment for testing
    data={
       "output_dir": "alignment_subworkflow_out",
       "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
       "run_id": 16906,
       "site_id": "CORI",
       "status": "uploading",
       "submission_id": "e7ef7456-bf05-4e12-ba3d-cd4e06727922",
       "upload_task_id": "444ac0b8-60f0-11eb-9905-0aa9ddbe2755"
    }
    """

    return data

def wait_for_run(env,run_id,check_tries,check_sleep):
    """ Wait for all the runs in run_ids list to finish."""
    tries = 1 
    while tries <= check_tries:
        # check whether the run has finished every 60 seconds
        cmd = "source ~/jaws-%s.sh > /dev/null && jaws run status %s" % (env,run_id)
        (r,o,e) = run(cmd)
        if r > 0:
            pytest.exit("stderr: %s" % e)

        status_output = json.loads(o)
        run_status = status_output["status"]

        if run_status == "download complete":
            return

        tries += 1
        time.sleep(check_sleep)

