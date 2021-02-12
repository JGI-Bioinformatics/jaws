import os
import json
import pytest
import smtplib
import time
import submission_utils as util


# @pytest.fixture(scope="module")
# def submit_wdl_and_wait():
#     """
#     This is a fixture that will submit a wdl for all functions to use.
#     This function returns the output of a wdl submission.
#     """
#     if os.path.exists(OUTDIR):
#         cmd = "rm -rf %s" % OUTDIR
#         (rc,stdout,stderr) = util.run(cmd)
#         if rc > 0:
#             os.exit("Failed to remove old output directory %s" % OUTDIR)
#
#     cmd = ". ~/jaws-%s.sh > /dev/null 2>&1 && jaws run submit %s %s %s %s" % (ENV,WDL,INPUT_JSON,OUTDIR,SITE)
#     (rc,stdout,stderr) = util.run(cmd)
#     if rc > 0:
#         pytest.exit("stderr: %s" % stderr)
#
#     assert rc == 0
#     data = json.loads(stdout)
#
#     # uncomment for testing
#     # """
#     # data={
#     #     "output_dir": "/global/cscratch1/sd/jfroula/JAWS/jaws/test/integration/end-to-end-tests/score-card-tests/out",
#     #     "output_endpoint": "9d6d994a-6d04-11e5-ba46-22000b92c6ec",
#     #     "run_id": 16679,
#     #     "site_id": "CORI",
#     #     "status": "uploading",
#     #     "submission_id": "6a16e3a3-a975-4924-b0e8-713215b8b771",
#     #     "upload_task_id": "444ac0b8-60f0-11eb-9905-0aa9ddbe2755"
#     # }
#     # """
#     run_id = str(data['run_id'])
#
#     # Wait for all the runs in run_ids list to finish.
#     check_tries=100
#     check_sleep=30
#     tries = 1
#     while tries <= check_tries:
#         # check whether the run has finished every 60 seconds
#         time.sleep(check_sleep)
#         cmd = "source ~/jaws-%s.sh > /dev/null && jaws run status %s" % (ENV,run_id)
#         (rc,stdout,stderr) = util.run(cmd)
#         if rc > 0:
#             pytest.exit("stderr: %s" % stderr)
#
#         status_output = json.loads(stdout)
#         run_status = status_output["status"]
#         result = status_output["result"]
#
#         if run_status == "download complete" and result == "succeeded":
#             return data
#
#         tries += 1
#
#     # if we got here the number of tries was exceeded, the run has not completed
#     pytest.exit("We have exeeded the wait time for the job to complete. You can increase the number of tries or sleep time.")

@pytest.fixture(scope="module")
def submit_fq_count_wdl(request):
    wdl = "fq_count.wdl"
    input_json = "fq_count.json"
    env = getattr(request.module, "env")
    outdir ="./fq_count_out"
    site = getattr(request.module, "site")
    util.submit_wdl(env, wdl, input_json, outdir, site)

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
        "--site",
        action="append",
        default=[],
        help="the JAWS site [cori|jgi] that will be used during submission",
    )

def pytest_generate_tests(metafunc):
    if "env" in metafunc.fixturenames:
        metafunc.parametrize("env", metafunc.config.getoption("env"))
    if "site" in metafunc.fixturenames:
        metafunc.parametrize("site", metafunc.config.getoption("site"))
