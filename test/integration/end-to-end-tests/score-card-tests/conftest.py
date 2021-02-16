import os
import json
import pytest
import smtplib
import time
import submission_utils

@pytest.fixture(scope="module")
def submit_fq_count_wdl(request):
    wdl = "./fq_count.wdl"
    input_json = "./fq_count.json"
    outdir ="./fq_count_out"

    # allow use to pass variables into the test functions via command line
    env = request.config.getoption("--env")
    site = request.config.getoption("--site")

    data = submission_utils.submit_wdl(env, wdl, input_json, outdir, site)
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
