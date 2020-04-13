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

    expected_amqp_sections = [
        ("host", "currenthost"),
        ("vhost", "current_vhost"),
        ("user", "daffy_duck"),
        ("password", "succotash"),
        ("queue", "rabbit_season"),
    ]

    expected_rpc_sections = [
        ("num_threads", "5"),
        ("max_retries", "3"),
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
        ("workflows_url", "http://localhost:8000/api/workflows/v1"),
        ("engine_status_url", "http://localhost:8000/engine/v1/status")]

    check_section("AMQP", expected_amqp_sections, cfg)
    check_section("RPC", expected_rpc_sections, cfg)
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

    check_section("AMQP", [("host", "https://rmq.nersc.gov")], cfg)
    check_section("AMQP", [("vhost", "/")], cfg)
    check_section("RPC", [("max_retries", "10")], cfg)
