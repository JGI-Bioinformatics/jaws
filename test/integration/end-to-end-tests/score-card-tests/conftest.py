import os
import json
import pytest
import smtplib
import time
import submission_utils as util

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
    outdir ="./fq_count_out"

    # allow use to pass variables into the test functions via command line
    env = request.config.getoption("--env")
    site = request.config.getoption("--site")

    data = util.submit_wdl(env, wdl, input_json, outdir, site)
    # print(data)  # used for debugging
    return data

@pytest.fixture(scope="session")
def submit_subworkflow_alignment(request):
    wdl = "./WDLs/jaws-alignment-example/main.wdl"
    input_json = "./WDLs/jaws-alignment-example/inputs.json"
    outdir ="./alignment_subworkflow_out"

    # allow use to pass variables into the test functions via command line
    env = request.config.getoption("--env")
    site = request.config.getoption("--site")

    data = util.submit_wdl(env, wdl, input_json, outdir, site)
    # print(data)  # used for debugging
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
