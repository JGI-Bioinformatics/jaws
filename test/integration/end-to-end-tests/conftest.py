import pytest
import submission_utils as util

check_tries = 50
check_sleep = 30


# @pytest.fixture(scope="session",autouse=True)
# def test_for_all_args(request):
#     target_dir = request.config.getoption("--dir")
#     site = request.config.getoption("--site")
#     if not target_dir or not site:
#         pytest.exit("Error: You are missing some arguments?\nUsage: \
#         pytest -n <number of tests in parallel> --capture=<[yes|no]> --verbose --dir <directory with tests> \
#         --site <[cori|jgi]> <directory or file>")


@pytest.fixture(scope="session")
def submit_fq_count_wdl(request):
    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")
    env = request.config.getoption("--env")

    wdl = target_dir + "/WDLs/fq_count.wdl"
    input_json = target_dir + "/test-inputs/fq_count.json"

    """
    data = {
    "output_dir": "/global/cfs/projectdirs/jaws/data-repository-staging/jfroula/CORI/6744",
    "run_id": 6744,
    "site_id": "CORI",
    "status": "uploading",
    "tag": ""
    }
    """

    data = util.submit_wdl(env, wdl, input_json, site)
    return data


@pytest.fixture(scope="session")
def submit_subworkflow_alignment(request):
    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")
    env = request.config.getoption("--env")

    wdl = target_dir + "/WDLs/jaws-alignment-example/main.wdl"
    input_json = target_dir + "/WDLs/jaws-alignment-example/inputs.json"

    data = util.submit_wdl(env, wdl, input_json, site)
    # data = {'output_dir': '/global/cfs/projectdirs/jaws/data-repository-staging/jfroula/JGI/6766', \
    # 'run_id': 6766, 'site_id': 'JGI', 'status': 'uploading'}

    # wait for run to complete
    id = data["run_id"]
    util.wait_for_run(id, env, check_tries, check_sleep)
    return data


@pytest.fixture(scope="session")
def submit_bad_task(request):
    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")
    env = request.config.getoption("--env")

    wdl = target_dir + "/WDLs/bad_task.wdl"
    input_json = target_dir + "/test-inputs/fq_count.json"

    data = util.submit_wdl_noexit(env, wdl, input_json, site)

    """
    data = {
    "output_dir": "/global/cfs/projectdirs/jaws/data-repository-staging/jfroula/CORI/6744",
    "run_id": 6744,
    "site_id": "CORI",
    "status": "uploading",
    "tag": ""
    }
    """

    id = data["run_id"]

    # wait for run to complete
    util.wait_for_run(id, env, check_tries, check_sleep)
    return data


@pytest.fixture(scope="session")
def submit_bad_docker(request):
    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")
    env = request.config.getoption("--env")

    wdl = target_dir + "/WDLs/bad_docker.wdl"
    input_json = target_dir + "/test-inputs/fq_count.json"

    data = util.submit_wdl(env, wdl, input_json, site)

    # wait for run to complete
    id = data["run_id"]
    util.wait_for_run(id, env, check_tries, check_sleep)

    # print(data)  # used for debugging
    return data


@pytest.fixture(scope="session")
def submit_skylake_250(request):

    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")
    env = request.config.getoption("--env")

    # skip this fixture if not run on cori
    if "cori" not in site.lower():
        pytest.skip("needs to run on cori only")

    wdl = target_dir + "/WDLs/skylake_test_250.wdl"
    input_json = target_dir + "/test-inputs/fq_count.json"
    site = "cori"

    data = util.submit_wdl(env, wdl, input_json, site)

    # wait for run to complete
    id = data["run_id"]
    util.wait_for_run(id, env, check_tries, check_sleep)

    return data


@pytest.fixture(scope="session")
def submit_skylake_500(request):
    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")
    env = request.config.getoption("--env")

    # skip this fixture if not run on cori
    if "cori" not in site.lower():
        pytest.skip("needs to run on cori only")

    wdl = target_dir + "/WDLs/skylake_test_500.wdl"
    input_json = target_dir + "/test-inputs/fq_count.json"

    data = util.submit_wdl(env, wdl, input_json, site)

    # wait for run to complete
    id = data["run_id"]
    util.wait_for_run(id, env, check_tries, check_sleep)

    return data


@pytest.fixture(scope="session")
def clone_tutorials_repo(request):
    # clone the jaws-tutorial-examples repo
    cmd = (
        "git clone "
        "https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git"
    )
    util.run(cmd)

    yield

    cmd = "rm -rf jaws-tutorial-examples/"
    util.run(cmd)

    yield

    cmd = "rm -rf jaws-tutorial-examples/"
    util.run(cmd)


# The addoption functions allows us to use flags to capture arguments on the command line.
def pytest_addoption(parser):
    parser.addoption(
        "--dir",
        action="store",
        default=["."],
        help="this is the path to where the WDLs and input.json files are.",
    )
    parser.addoption(
        "--site",
        action="store",
        help="the JAWS site [cori|jgi] that will be used during submission",
    )
    parser.addoption(
        "--env",
        action="store",
        help="the JAWS environment [dev|staging|prod] that will be used during submission",
    )


# These functions allows an argument to be passed into the test functions
@pytest.fixture
def dir(request):
    return request.config.getoption("--dir")


@pytest.fixture
def site(request):
    return request.config.getoption("--site")


@pytest.fixture
def env(request):
    return request.config.getoption("--env")
