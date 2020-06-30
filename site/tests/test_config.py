import jaws_site.config


def check_section(section_to_test, expected_entries, actual_config):
    for section, expected in expected_entries:
        assert actual_config.get(section_to_test, section) == expected


def test_overwrite_all_default_values(config_file):

    # Because Configuration is a singleton, we call a destructor method to
    # remove any old reference.
    jaws_site.config.Configuration._destructor()

    config_path = config_file
    cfg = jaws_site.config.Configuration(config_path)

    expected_jtm_rpc_server_sections = [
        ("host", "localhost"),
        ("vhost", "site_vhost"),
        ("user", "jaws"),
        ("password", "passw0rd1"),
        ("num_threads", "5"),
        ("max_retries", "3"),
    ]

    expected_central_rpc_server_sections = [
        ("host", "currenthost"),
        ("vhost", "eagle"),
        ("user", "daffy_duck"),
        ("password", "succotash"),
        ("num_threads", "5"),
        ("max_retries", "3"),
    ]

    expected_rpc_client_sections = [
        ("host", "currenthost"),
        ("vhost", "central_vhost"),
        ("user", "daffy_duck"),
        ("password", "succotash"),
    ]

    expected_globus_sections = [
        ("client_id", "foghorn_leghorn"),
        ("endpoint_id", "rooster"),
        ("root_dir", "cwd")
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
        ("staging_subdirectory", "staging"),
        ("results_subdirectory", "results"),
    ]

    expected_cromwell_sections = [
        ("url", "http://localhost:8000")]

    check_section("JTM_RPC_SERVER", expected_jtm_rpc_server_sections, cfg)
    check_section("CENTRAL_RPC_SERVER", expected_central_rpc_server_sections, cfg)
    check_section("CENTRAL_RPC_CLIENT", expected_rpc_client_sections, cfg)
    check_section("GLOBUS", expected_globus_sections, cfg)
    check_section("DB", expected_db_sections, cfg)
    check_section("SITE", expected_site_sections, cfg)
    check_section("CROMWELL", expected_cromwell_sections, cfg)


def test_config_overwrite_partial_values(partial_config):

    # Because Configuration is a singleton, we call a destructor method to
    # remove any old reference.
    jaws_site.config.Configuration._destructor()

    config_path = partial_config
    cfg = jaws_site.config.Configuration(config_path)

    check_section("CENTRAL_RPC_SERVER", [("host", "https://rmq.nersc.gov")], cfg)
    check_section("CENTRAL_RPC_SERVER", [("vhost", "site")], cfg)
    check_section("CENTRAL_RPC_SERVER", [("max_retries", "10")], cfg)
