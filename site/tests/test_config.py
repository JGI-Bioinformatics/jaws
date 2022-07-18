import jaws_site.config
import os


# Move large variable declarations out method definition to here
expected_local_rpc_server_sections = [
    ("host", "localhost"),
    ("vhost", "jaws_test"),
    ("queue", "site_rpc"),
    ("user", "jaws"),
    ("password", "passw0rd1"),
    ("password2", '${LOCAL_RPC_SERVER_PASSWORD}'),
    ("password3", '${LOCAL_RPC_SERVER_PASSWORD}$'),
    ("num_threads", "5"),
    ("max_retries", "3"),
]

expected_local_rpc_server_sections_env = [
    ("host", "localhost"),
    ("vhost", "jaws_test"),
    ("queue", "site_rpc"),
    ("user", "jaws"),
    ("password", "passw0rd1"),
    ("password2", 'password'),
    ("password3", 'password$'),
    ("num_threads", "5"),
    ("max_retries", "3"),
]

expected_central_rpc_server_sections = [
    ("host", "currenthost"),
    ("vhost", "jaws_test"),
    ("user", "jaws_eagle"),
    ("password", "succotash"),
    ("num_threads", "5"),
    ("max_retries", "3"),
]

expected_rpc_client_sections = [
    ("host", "currenthost"),
    ("vhost", "jaws_test"),
    ("queue", "central_rpc"),
    ("user", "jaws_eagle"),
    ("password", "succotash"),
]

expected_db_sections = [
    ("dialect", "mysql+mysqlconnector"),
    ("host", "myhost"),
    ("port", "60032"),
    ("user", "elmer_fudd"),
    ("password", "hunting"),
    ("db", "hunting_sites"),
]

expected_site_sections = [
    ("id", "eagle"),
    ("inputs_dir", "/global/scratch/jaws/jaws-dev/inputs"),
]

expected_cromwell_sections = [
    ("url", "http://localhost:8000")]


def check_section(section_to_test, expected_entries, actual_config):
    for section, expected in expected_entries:
        print(section, expected, actual_config.get(section_to_test, section))
        assert actual_config.get(section_to_test, section) == expected


def test_overwrite_all_default_values(config_file):

    # Because Configuration is a singleton, we call a destructor method to
    # remove any old reference.
    jaws_site.config.Configuration._destructor()

    config_path = config_file
    cfg = jaws_site.config.Configuration(config_path)

    print("Checking LOCAL_RPC_SERVER without environment variables set.")
    # Clear potentially conflicting env vars
    try:
        del os.environ['LOCAL_RPC_SERVER_PASSWORD']
    except KeyError:
        pass
    check_section("LOCAL_RPC_SERVER", expected_local_rpc_server_sections, cfg)
    print("Checking LOCAL_RPC_SERVER with environment variables set.")
    os.environ['LOCAL_RPC_SERVER_PASSWORD'] = 'password'
    check_section("LOCAL_RPC_SERVER", expected_local_rpc_server_sections_env, cfg)
    check_section("RPC_SERVER", expected_central_rpc_server_sections, cfg)
    check_section("CENTRAL_RPC_CLIENT", expected_rpc_client_sections, cfg)
    # Clear potentially conflicting env vars
    try:
        del os.environ['PROJECT_NAME']
    except KeyError:
        pass
    os.environ['PROJECT_NAME'] = 'jaws'
    check_section("DB", expected_db_sections, cfg)
    check_section("SITE", expected_site_sections, cfg)
    check_section("CROMWELL", expected_cromwell_sections, cfg)


def test_config_overwrite_partial_values(partial_config):

    # Because Configuration is a singleton, we call a destructor method to
    # remove any old reference.
    jaws_site.config.Configuration._destructor()

    config_path = partial_config
    cfg = jaws_site.config.Configuration(config_path)

    check_section("RPC_SERVER", [("host", "https://rmq.nersc.gov")], cfg)
    check_section("RPC_SERVER", [("vhost", "jaws_test")], cfg)
    check_section("RPC_SERVER", [("max_retries", "10")], cfg)
