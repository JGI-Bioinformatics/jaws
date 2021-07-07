import sys
import json
import pytest
import time
import functools
from subprocess import Popen, PIPE

# flush the printstream so that the output can be redirected to a file to prevent
# loss of test results if the ssh connection is lost
print = functools.partial(print, flush=True)


def run(cmd):
    output = Popen(cmd, stdout=PIPE,
                   stderr=PIPE, shell=True,
                   universal_newlines=True)

    stdout, stderr = output.communicate()
    rc = output.returncode
    if rc:
        sys.stderr.write("Error: This command returned an error code greater than 0: \n%s\n%s" % (cmd, stderr))

    return rc, stdout, stderr


def submit_wdl(env, wdl, input_json, site, tag="submit_fq_count_wdl"):
    """
    This is a fixture that will submit a wdl for all functions to use.
    This function returns the output of a wdl submission.
    """
    # the pipe > /dev/null 2>&1 is needed below because otherwise the info printed from the
    # activation command causes an error when we try to do json load later

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws submit --tag %s --no-cache %s %s %s" % (
        env, tag, wdl, input_json, site)
    (rc, stdout, stderr) = run(cmd)
    if rc > 0:
        if stderr:
            pytest.exit("stderr: %s" % stderr)
        if stdout:
            pytest.exit("stdout: %s" % stdout)
        else:
            pytest.exit(
                "The return code from the command was greater than 0. There was no stderr accompanying this failure. "
                "\n%s" % cmd)

    assert rc == 0
    data = json.loads(stdout)
    return data


def submit_wdl_noexit(env, wdl, input_json, site):
    """
    This is a fixture that will submit a wdl that is expected to error out. 
    I will not exit the function if there is an error, but will just return the stderr, stdout, and rc.
    """

    # the pipe > /dev/null 2>&1 is needed below because otherwise the info printed from the
    # activation command causes an error when we try to do json load later

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws submit --tag submit_wdl_noexit --no-cache %s %s %s" % (
        env, wdl, input_json, site)
    (rc, stdout, stderr) = run(cmd)
    data = json.loads(stdout)

    return data


def wait_for_run(run_id, env, check_tries, check_sleep):
    """ Wait for the run to finish."""
    tries = 1
    run_status = "unknown"
    while tries <= check_tries:
        previous_status = run_status

        # check whether the run has finished
        cmd = "source ~/jaws-%s.sh > /dev/null && jaws status %s" % (env, run_id)
        (r, o, e) = run(cmd)
        if r > 0:
            pytest.exit("stderr: %s" % e)

        status_output = json.loads(o)
        run_status = status_output["status"]

        # if the run has finished we are done waiting
        if run_status == "download complete":
            return

        # if the run status has changed, reset the tries counter to 1
        if run_status != previous_status:
            tries = 1
        # otherwise the run status is the same as last time, so increment the runs counter
        else:
            tries += 1

        # sleep for awhile before checking the status again
        time.sleep(check_sleep)

    # if we got here then the number of tries has been exceeded, but the run is still not finished
    error_message = "The test has timed out while waiting for run %s to complete" % run_id
    raise Exception(error_message)


def wait_for_run_and_check_for_success(run_id, env, check_tries, check_sleep):
    wait_for_run(run_id, env, check_tries, check_sleep)

    cmd = "source ~/jaws-%s.sh > /dev/null && jaws status %s" % (env, run_id)
    (rc, stdout, stderr) = run(cmd)

    status_info = json.loads(stdout)
    assert status_info["status"] == "download complete", \
        "\n** Run %s took too long - last state seen: %s" % (run_id, status_info["status"])
    assert status_info["result"] == "succeeded", \
        "\n** Run %s did not succeed - status was: %s" % (run_id, status_info["result"])


def timestamp_dir(base_dir):
    # create timestamp string to make output directory unique
    ts = str(time.time()).replace(".", "")
    if base_dir[-1] != '/':
        base_dir = base_dir + '/'
    return base_dir + ts
