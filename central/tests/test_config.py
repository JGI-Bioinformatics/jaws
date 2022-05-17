import jaws_central.config
import os


def check_section(section_to_test, expected_entries, actual_config):
    for key, expected in expected_entries:
        assert actual_config.get(section_to_test, key) == expected


def check_site(site_to_test, expected_entries, actual_config):
    for key, expected in expected_entries:
        assert actual_config.get_site_param(site_to_test, key) == expected


def check_site_info(site_to_test, expected_entries, actual_config):
    site_params = actual_config.get_site_info(site_to_test)
    for key, expected in expected_entries:
        assert key in site_params
        assert site_params[key] == expected


def test_check_all_values(config_file):

    config_path = config_file
    cfg = jaws_central.config.Configuration(config_path)
    # Make sure that the _vars attribute in the interpolation object is "NONE" because
    # there was no env_override value set
    print("Verifying that the _vars instance variable is NONE")
    assert cfg.config._interpolation._vars is None

    expected_db_parameters = [
        ("dialect", "mysql+mysqlconnector"),
        ("host", "db.foo.com"),
        ("port", "3306"),
        ("user", "jaws"),
        ("password", "passw0rd1"),
        ("db", "jaws"),
        ("host2", "${JAWS_DB_HOST}"),
        ("password2", "${JAWS_DB_PASSWORD}123"),
    ]

    expected_db_parameters_env = [
        ("dialect", "mysql+mysqlconnector"),
        ("host", "db.foo.com"),
        ("port", "3306"),
        ("user", "jaws"),
        ("password", "passw0rd1"),
        ("db", "jaws"),
        ("host2", "db.foobar.com"),
        ("password2", "password123"),
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
        ("inputs_dir", "/global/scratch/jaws/jaws-dev/inputs"),
        ("inputs_dir2", "${SCRATCH_ROOT}/jaws/jaws-dev/inputs"),
        ("max_ram_gb", "1024"),
    ]

    expected_site2_parameters = [
        ("site_id", "NERSC"),
        ("globus_endpoint", "YYYY"),
        ("globus_host_path", "/"),
        ("inputs_dir", "/global/cscratch/sd1/jaws/jaws-dev/inputs"),
        ("max_ram_gb", "2048"),
    ]

    print("Checking DB section without environment variables set")
    # Clear any conflicting keys in the environment
    try:
        del os.environ["JAWS_DB_HOST"]
        del os.environ["JAWS_DB_PASSWORD"]
    except KeyError:
        pass
    check_section("DB", expected_db_parameters, cfg)
    print("Checking DB section with environment variables set")
    os.environ["JAWS_DB_HOST"] = "db.foobar.com"
    os.environ["JAWS_DB_PASSWORD"] = "password"
    check_section("DB", expected_db_parameters_env, cfg)
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
