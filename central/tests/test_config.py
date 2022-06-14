import jaws_central.config
import os
import pytest


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
    try:
        jaws_central.config.Configuration._destructor()
    except Exception:
        pass
    config_path = config_file
    print("Test exception when env_override value is longer than 20 chars")
    with pytest.raises(ValueError):
        jaws_central.config.Configuration(config_path, "012345678901234567890_")

    try:
        jaws_central.config.Configuration._destructor()
    except Exception:
        pass
    print("Test exception when env_override value is shorter than 3 chars")
    with pytest.raises(ValueError):
        jaws_central.config.Configuration(config_path, "0_")

    try:
        jaws_central.config.Configuration._destructor()
    except Exception:
        pass

    cfg = jaws_central.config.Configuration(config_path)
    # Make sure that the _vars attribute in the interpolation object is "NONE" because
    # there was no env_override value set
    print("Verifying that the _vars instance variable is NONE")
    assert cfg.config._vars is None
    print("Instances of jaws_central.config.Configuration: ", dict(jaws_central.config.Singleton._instances))
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

    expected_db_parameters_env_override = [
        ("dialect", "mysql+mysqlconnector"),
        ("host", "db.foo.com"),
        ("port", "3306"),
        ("user", "jaws"),
        ("password", "passw0rd1"),
        ("db", "jaws"),
        ("host2", "db_over.foobar.com"),
        ("password2", "over_password123"),
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

    expected_site1_parameters2 = [
        ("host", "jgi-host.foobar.org"),
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

    expected_site2_parameters2 = [
        ("host", "nersc-host.foobar.org"),
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

    try:
        jaws_central.config.Configuration._destructor()
    except Exception:
        pass

    os.environ["ENV__host2"] = "db_over.foobar.com"
    os.environ["ENV__password2"] = "over_password123"
    print("Recreating new JAWSConfig object with env_override set to ENV__")
    cfg = jaws_central.config.Configuration(config_path, "ENV__")
    print("Instances of jaws_central.config.Configuration: ", dict(jaws_central.config.Singleton._instances))
    print("Verifying that _vars attribute has picked up environment variables")
    assert(cfg.config._vars is not None)
    assert(len(cfg.config._vars) == 2)
    assert(cfg.config._vars['host2'] == "db_over.foobar.com")
    assert(cfg.config._vars['password2'] == "over_password123")

    print("Retesting with env_prefix set")
    check_section("DB", expected_db_parameters_env_override, cfg)
    check_section("GLOBUS", expected_globus_parameters, cfg)
    check_site("JGI", expected_site1_parameters, cfg)
    check_site_info("NERSC", expected_site2_parameters, cfg)

    try:
        jaws_central.config.Configuration._destructor()
    except Exception:
        pass
    print("Retesting with env_prefix set and section names")
    os.environ["ENV__SITE:JGI_host"] = "jgi-host.foobar.org"
    os.environ["ENV__SITE:NERSC_host"] = "nersc-host.foobar.org"

    # Verify we get a KeyError when mixing section names and non-section
    with pytest.raises(KeyError):
        jaws_central.config.Configuration(config_path, "ENV__")
    del(os.environ["ENV__host2"])
    del(os.environ["ENV__password2"])

    try:
        jaws_central.config.Configuration._destructor()
    except Exception:
        pass
    # We've cleared the key error, see if we've SITE:JGI_host and
    # SITE:NERSC_host settings are in place
    cfg = jaws_central.config.Configuration(config_path, "ENV__")
    assert("SITE:JGI" in cfg.config._vars)
    assert("SITE:NERSC" in cfg.config._vars)
    assert(cfg.config['SITE:JGI']['host'] == os.environ["ENV__SITE:JGI_host"])
    assert(cfg.config['SITE:NERSC']['host'] == os.environ["ENV__SITE:NERSC_host"])
    #  check_section("GLOBUS", expected_globus_parameters, cfg)
    # check_site("JGI", expected_site1_parameters2, cfg)
    # check_site_info("NERSC", expected_site2_parameters2, cfg)


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
