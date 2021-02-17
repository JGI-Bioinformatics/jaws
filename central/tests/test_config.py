import jaws_central.config


def check_section(section_to_test, expected_entries, actual_config):
    for key, expected in expected_entries:
        assert actual_config.get(section_to_test, key) == expected


def check_site(site_to_test, expected_entries, actual_config):
    for key, expected in expected_entries:
        assert actual_config.get_site(site_to_test, key) == expected


def check_site_info(site_to_test, expected_entries, actual_config):
    site_params = actual_config.get_site_info(site_to_test)
    for key, expected in expected_entries:
        assert key in site_params
        assert site_params[key] == expected


def test_check_all_values(config_file):

    config_path = config_file
    cfg = jaws_central.config.Configuration(config_path)

    expected_db_parameters = [
        ("dialect", "mysql+mysqlconnector"),
        ("host", "db.foo.com"),
        ("port", "3306"),
        ("user", "jaws"),
        ("password", "passw0rd1"),
        ("db", "jaws"),
    ]

    expected_globus_parameters = [
        ("client_id", "ZZZZ"),
        ("client_secret", "AAAAA"),
    ]

    expected_site1_parameters = [
        ("host", "rmq.jaws.gov"),
        ("user", "jaws"),
        ("password", "passw0rd3"),
        ("vhost", "jaws"),
        ("globus_endpoint", "XXXX"),
        ("globus_host_path", "/global/scratch/jaws"),
        ("uploads_subdir", "jaws-dev/uploads"),
        ("max_ram_gb", "1024"),
    ]

    expected_site2_parameters = [
        ("site_id", "NERSC"),
        ("globus_endpoint", "YYYY"),
        ("globus_host_path", "/"),
        ("uploads_subdir", "/global/cscratch/sd1/jaws/jaws-dev/uploads"),
        ("max_ram_gb", "2048"),
    ]

    check_section("DB", expected_db_parameters, cfg)
    check_section("GLOBUS", expected_globus_parameters, cfg)
    check_site("JGI", expected_site1_parameters, cfg)
    check_site_info("NERSC", expected_site2_parameters, cfg)


def test_config_overwrite_partial_values(partial_config):

    # call destructor to remove old references
    jaws_central.config.Configuration._destructor()

    config_path = partial_config
    cfg = jaws_central.config.Configuration(config_path)

    expected_db_sections = [
        ("dialect", "mysql+mysqlconnector"),
        ("host", "db.foo.com"),
        ("port", "3305"),
        ("user", "j4ws"),
        ("password", "p455w0rd1"),
        ("db", "jaws"),
    ]

    expected_globus_sections = [
        ("client_id", "AAAA"),
    ]

    check_section("DB", expected_db_sections, cfg)
    check_section("GLOBUS", expected_globus_sections, cfg)
