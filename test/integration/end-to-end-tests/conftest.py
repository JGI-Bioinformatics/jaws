import pytest
import submission_utils as util

check_sleep = 60
check_tries = 360


@pytest.fixture(scope="session")
def submit_fq_count_wdl(request):
    """ This will submit fq_count.wdl and fq_count.json and NOT wait for it to complete.

    :param a request object for capturing CLI arguments
    :type object 
    :rtype dictionary 
    :return output from jaws submit 
    """

    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")

    wdl = target_dir + "/WDLs/fq_count.wdl"
    input_json = target_dir + "/test-inputs/fq_count.json"
    """
    data =  {
    "max_ram_gb": 10,
    "run_id": 7105,
    }
    """
    data = util.submit_wdl(wdl, input_json, site)
    return data


@pytest.fixture(scope="session")
def submit_subworkflow_alignment(request):
    """ This will submit a simple WDL that includes a subworkflow.

    :param a request object for capturing CLI arguments
    :type object 
    :rtype dictionary 
    :return output from jaws submit 
    """
    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")

    wdl = target_dir + "/WDLs/jaws-alignment-example/main.wdl"
    input_json = target_dir + "/WDLs/jaws-alignment-example/inputs.json"

    data = util.submit_wdl(wdl, input_json, site)

    # wait for run to complete
    id = data["run_id"]
    util.wait_for_run(id, check_tries, check_sleep)
    return data


@pytest.fixture(scope="session")
def submit_bad_task(request):
    """ This will submit a fq_count.wdl WDL that has a bad bash command in the command{} section. This is for testing the error logs.

    :param a request object for capturing CLI arguments
    :type object 
    :rtype dictionary 
    :return output from jaws submit 
    """
    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")

    wdl = target_dir + "/WDLs/bad_task.wdl"
    input_json = target_dir + "/test-inputs/fq_count.json"

    data = util.submit_wdl_noexit(wdl, input_json, site)

    id = data["run_id"]

    # wait for run to complete
    util.wait_for_run(id, check_tries, check_sleep=30)
    return data


@pytest.fixture(scope="session")
def submit_bad_docker(request):
    """ This will submit a fq_count.wdl WDL that has a non-existent docker image name in the runtime{} section. This is for testing the error logs.

    :param a request object for capturing CLI arguments
    :type object 
    :rtype dictionary 
    :return output from jaws submit 
    """
    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")

    wdl = target_dir + "/WDLs/bad_docker.wdl"
    input_json = target_dir + "/test-inputs/fq_count.json"

    data = util.submit_wdl( wdl, input_json, site)

    # wait for run to complete
    id = data["run_id"]
    util.wait_for_run(id, check_tries, check_sleep)

    # print(data)  # used for debugging
    return data


@pytest.fixture(scope="session")
def submit_bad_sub_task(request):
    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")

    wdl = target_dir + "/WDLs/main_bad_sub_task.wdl"
    input_json = target_dir + "/test-inputs/main_bad_sub_task.json"

    data = util.submit_wdl_noexit(wdl, input_json, site)

    id = data["run_id"]

    # wait for run to complete
    util.wait_for_run(id, check_tries, check_sleep=30)
    return data


@pytest.fixture(scope="session")
def submit_scatter_timeout(request):
    # allow user to pass variables into the test functions via command line
    target_dir = request.config.getoption("--dir")
    site = request.config.getoption("--site")

    wdl = target_dir + "/../../../../examples/leo_dapseq/leo_15_min.wdl"
    input_json = target_dir + "/../../../../examples/leo_dapseq/shortened-100.json"

    data = util.submit_wdl_noexit(wdl, input_json, site)

    id = data["run_id"]

    # wait for run to complete
    util.wait_for_run(id, check_tries, check_sleep)
    return data


@pytest.fixture(scope="session")
def clone_tutorials_repo(request):
    """ Clones https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples.

    :param a request object for capturing CLI arguments
    :type object 
    :rtype none 
    :return none
    """
    # clone the jaws-tutorial-examples repo
    cmd = (
        "git clone https://code.jgi.doe.gov/official-jgi-workflows/jaws-tutorial-examples.git"
    )
    util.run(cmd)

    yield

    cmd = "rm -rf jaws-tutorial-examples/"
    util.run(cmd)


# The addoption functions allows us to use flags to capture arguments on the command line.
def pytest_addoption(parser):
    """ The parser.addoption function allows us to use flags to capture CLI arguments that can then be used in our test functions as if they were fixtures. 

    :param a request object for capturing CLI arguments
    :type object
    :rtype none
    :return none
    """
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

# These functions allows an argument to be passed into the test functions
@pytest.fixture
def dir(request):
    """ This CLI argument is for the parent directory in which the WDL & input.json files live."""
    return request.config.getoption("--dir")


@pytest.fixture
def site(request):
    """ The --site CLI argument is for passing the site [cori|jgi] to test the functions"""
    return request.config.getoption("--site")

