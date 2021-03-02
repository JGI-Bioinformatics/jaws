import os
import json
import pytest
import smtplib
import time
import submission_utils as util

check_tries=50
check_sleep=30

@pytest.fixture(scope="session",autouse=True)
def test_for_all_args(request):
    env = request.config.getoption("--env")
    site = request.config.getoption("--site")
    if not env or not site: 
        pytest.exit("Error: You are missing some arguments?\nUsage: pytest -n <number of tests in parallel> --capture=<[yes|no]> --verbose --env <[prod|staging|dev]> --site <[cori|jgi]> <directory or file>")

@pytest.fixture(scope="session")
def submit_fq_count_wdl(request):
    wdl = "./WDLs/fq_count.wdl"
    input_json = "./test-inputs/fq_count.json"

    # allow user to pass variables into the test functions via command line
    env = request.config.getoption("--env")
    site = request.config.getoption("--site")

    data = util.submit_wdl(env, wdl, input_json, site)
    return data

@pytest.fixture(scope="session")
def submit_subworkflow_alignment(request):
    wdl = "./WDLs/jaws-alignment-example/main.wdl"
    input_json = "./WDLs/jaws-alignment-example/inputs.json"

    # allow user to pass variables into the test functions via command line
    env = request.config.getoption("--env")
    site = request.config.getoption("--site")

    data = util.submit_wdl(env, wdl, input_json, site)

    # wait for run to complete
    run_id = data['run_id']
    util.wait_for_run(env,run_id,check_tries,check_sleep)
    # print(data)  # used for debugging
    return data


@pytest.fixture(scope="session")
def submit_bad_task(request):
    wdl = "./WDLs/bad_task.wdl"
    input_json = "./test-inputs/fq_count.json"

    # allow user to pass variables into the test functions via command line
    env = request.config.getoption("--env")
    site = request.config.getoption("--site")

    data = util.submit_wdl_noexit(env, wdl, input_json, site)
    run_id = data['run_id']

    # wait for run to complete
    util.wait_for_run(env,run_id,check_tries,check_sleep)
    return data


@pytest.fixture(scope="session")
def submit_bad_docker(request):
    wdl = "./WDLs/bad_docker.wdl"
    input_json = "./test-inputs/fq_count.json"

    # allow user to pass variables into the test functions via command line
    env = request.config.getoption("--env")
    site = request.config.getoption("--site")

    data = util.submit_wdl(env, wdl, input_json, site)

    # wait for run to complete
    run_id = data['run_id']
    util.wait_for_run(env,run_id,check_tries,check_sleep)

    # print(data)  # used for debugging
    return data

@pytest.fixture(scope="session")
def submit_skylake_250(request):
    wdl = "./WDLs/skylake_test_250.wdl"
    input_json = "./test-inputs/fq_count.json"
    site="cori"

    # allow user to pass variables into the test functions via command line
    env = request.config.getoption("--env")

    data = util.submit_wdl(env, wdl, input_json, site)

    # wait for run to complete
    run_id = data['run_id']
    util.wait_for_run(env,run_id,check_tries,check_sleep)

    return data

@pytest.fixture(scope="session")
def submit_skylake_500(request):
    wdl = "./WDLs/skylake_test_500.wdl"
    input_json = "./test-inputs/fq_count.json"
    site="cori"

    # allow user to pass variables into the test functions via command line
    env = request.config.getoption("--env")

    data = util.submit_wdl(env, wdl, input_json, site)

    # wait for run to complete
    run_id = data['run_id']
    util.wait_for_run(env,run_id,check_tries,check_sleep)

    return data


#
# The next two functions allows us to use the --env to capture the environment [prod|staging|dev]. 
# This environment is an argument that can be passed into the test functions
#
def pytest_addoption(parser):
    parser.addoption(
        "--env",
        action="store",
        default=[],
        help="testing environment [prod|staging|dev] passed to test functions",
    )
    parser.addoption(
        "--site",
        action="store",
        default=[],
        help="the JAWS site [cori|jgi] that will be used during submission",
    )

@pytest.fixture
def env(request):
    return request.config.getoption("--env")

@pytest.fixture
def site(request):
    return request.config.getoption("--site")


# def pytest_generate_tests(metafunc):
#     if "env" in metafunc.fixturenames:
#         metafunc.parametrize("env", metafunc.config.getoption("env"))
#     if "site" in metafunc.fixturenames:
#         metafunc.parametrize("site", metafunc.config.getoption("site"))
