import os
import json
import pytest
import smtplib
import time
from subprocess import Popen, PIPE
import configparser

# set some environmental vars
config = configparser.ConfigParser()
config.read(os.environ.get('MYINI_FILE'))

WDL        = config['wdl']['wdl']
INPUT_JSON = config['wdl']['input_json']
ENV        = config['wdl']['env']
OUTDIR     = config['wdl']['outdir']
SITE       = config['wdl']['site']


def run(cmd):
    output = Popen(cmd, stdout=PIPE,
             stderr=PIPE, shell=True,
             universal_newlines=True)

    stdout,stderr=output.communicate()
    rc=output.returncode

    return rc,stdout,stderr

@pytest.fixture(scope="module")
def submit_wdl_and_wait():
    """
    This is a fixture that will submit a wdl for all functions to use.  
    This function returns the output of a wdl submission. 
    """
    cmd = ". ~/jaws-%s.sh > /dev/null 2>&1 && jaws run submit %s %s fq_count_out cori" % (env,wdl,inputs)
    (rc,stdout,stderr) = run(cmd)
    if rc > 0:
        pytest.exit("stderr: %s" % stderr)

    assert rc == 0
    data = json.loads(stdout)

    # uncomment for testing
    """
    data={
        "output_dir": "/global/cscratch1/sd/jfroula/JAWS/jaws/test/integration/end-to-end-tests/score-card-tests/out",
        "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
        "run_id": 16679,
        "site_id": "CORI",
        "status": "uploading",
        "submission_id": "6a16e3a3-a975-4924-b0e8-713215b8b771",
        "upload_task_id": "444ac0b8-60f0-11eb-9905-0aa9ddbe2755"
    }
    """
    run_id = str(data['run_id'])

    # Wait for all the runs in run_ids list to finish.
    check_tries=100
    check_sleep=30
    tries = 1
    while tries <= check_tries:
        # check whether the run has finished every 60 seconds
        time.sleep(check_sleep)
        cmd = "source ~/jaws-%s.sh > /dev/null && jaws run status %s" % (ENV,run_id)
        (rc,stdout,stderr) = run(cmd)
        if rc > 0:
            pytest.exit("stderr: %s" % stderr)

        status_output = json.loads(stdout)
        run_status = status_output["status"]
        result = status_output["result"]

        if run_status == "download complete" and result == "succeeded":
            return data

        tries += 1

    # if we got here the number of tries was exceeded, the run has not completed
    pytest.exit("tries exceeded")

@pytest.fixture(scope="module")
def submit_wdl():
    """
    This is a fixture that will submit a wdl for all functions to use.  
    This function returns the output of a wdl submission. 
    """
    if os.path.exists(OUTDIR):
        cmd = "rm -rf %s" % OUTDIR
        (rc,stdout,stderr) = run(cmd)
        if rc > 0:
            os.exit("Failed to remove old output directory %s" % OUTDIR)

    cmd = ". ~jfroula/jaws-%s.sh > /dev/null 2>&1 && jaws run submit %s %s fq_count_out cori" % (ENV,WDL,INPUT_JSON)
    (rc,stdout,stderr) = run(cmd)
    if rc > 0:
        pytest.exit("stderr: %s" % stderr)
    
    assert rc == 0
    data = json.loads(stdout)

    # uncomment for testing
    #data={
    #    "output_dir": "fq_count_out",
    #    "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
    #    "run_id": 16405,
    #    "site_id": "CORI",
    #    "status": "uploading",
    #    "submission_id": "47a555c7-07a6-442c-a2f1-d0319f2e3008",
    #    "upload_task_id": "444ac0b8-60f0-11eb-9905-0aa9ddbe2755"
    #}
    return data
#
# The next two functions allows us to use the --env to capture the environment [prod|staging|dev]. 
# This environment is an argument that can be passed into the test functions
#
def pytest_addoption(parser):
    parser.addoption(
        "--env",
        action="append",
        default=[],
        help="testing environment [prod|staging|dev] passed to test functions",
    )
    parser.addoption(
        "--wdl",
        action="append",
        default=[],
        help="the wdl that will be submitted",
    )

def pytest_generate_tests(metafunc):
    if "env" in metafunc.fixturenames:
        metafunc.parametrize("env", metafunc.config.getoption("env"))
    if "wdl" in metafunc.fixturenames:
        metafunc.parametrize("wdl", metafunc.config.getoption("wdl"))
