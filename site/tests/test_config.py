import os

import jaws_site.config
import pytest
from jaws_site.config import ConfigurationError

# Move large variable declarations out method definition to here
expected_rmq_sections = [
    ("host", "localhost"),
    ("vhost", "jaws_test"),
    ("user", "jaws"),
    ("password", "password"),
    ("num_threads", "5"),
    ("max_retries", "3"),
]

expected_db_sections = [
    ("dialect", "mysql+mysqlconnector"),
    ("host", "localhost"),
    ("port", "3306"),
    ("user", "elmer_fudd"),
    ("password", "hunting"),
    ("db", "hunting_sites"),
]

expected_db_sections_env_override = [
    ("dialect", "mysql+mysqlconnector"),
    ("host", "localhost"),
    ("port", "3306"),
    ("user", "elmer_fudd"),
    ("password", "hunting"),
    ("db", "hunting_sites"),
    ("host2", "db_over.foobar.com"),
    ("password2", "over_password123"),
]

expected_db_sections_env_override2 = [
    ("dialect", "mysql+mysqlconnector"),
    ("host", "db-host.foo.com"),
    ("port", "3306"),
    ("user", "elmer_fudd"),
    ("password", "hunting"),
    ("db", "hunting_sites"),
    ("host2", "db.foobar.com"),
    ("password2", "password123"),
]

expected_site_sections = [
    ("id", "eagle"),
    ("inputs_dir", "/global/scratch/jaws/jaws-dev/inputs"),
]

expected_cromwell_sections = [("url", "http://localhost:8000")]


def check_section(section_to_test, expected_entries, actual_config):
    for section, expected in expected_entries:
        # print(section, expected, actual_config.get(section_to_test, section))
        assert actual_config.get(section_to_test, section) == expected


def test_overwrite_all_default_values(config_file):
    # Because Configuration is a singleton, we call a destructor method to
    # remove any old reference.
    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass
    config_path = config_file
    cfg = jaws_site.config.Configuration(config_path)

    # print("Checking RMQ without environment variables set.")
    # Clear potentially conflicting env vars
    try:
        del os.environ["RMQ_PASSWORD"]
    except KeyError:
        pass
    # print("Checking RMQ with environment variables set.")
    os.environ["RMQ_PASSWORD"] = "password"
    check_section("RMQ", expected_rmq_sections, cfg)
    # Clear potentially conflicting env vars
    try:
        del os.environ["PROJECT_NAME"]
    except KeyError:
        pass
    os.environ["PROJECT_NAME"] = "jaws"
    check_section("DB", expected_db_sections, cfg)
    check_section("SITE", expected_site_sections, cfg)
    check_section("CROMWELL", expected_cromwell_sections, cfg)


def test_config_overwrite_partial_values(partial_config):
    # Because Configuration is a singleton, we call a destructor method to
    # remove any old reference.
    jaws_site.config.Configuration._destructor()

    config_path = partial_config
    cfg = jaws_site.config.Configuration(config_path)

    check_section("RMQ", [("host", "https://rmq.nersc.gov")], cfg)
    check_section("RMQ", [("vhost", "jaws_test")], cfg)
    check_section("RMQ", [("max_retries", "3")], cfg)


def test_env_override(config_file):
    """
    Separate checks for the environment override code
    """
    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass
    config_path = config_file
    # print("Test exception when env_override value is longer than 20 chars")
    with pytest.raises(ValueError):
        jaws_site.config.Configuration(config_path, "012345678901234567890_")

    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass
    # print("Test exception when env_override value is shorter than 3 chars")
    with pytest.raises(ValueError):
        jaws_site.config.Configuration(config_path, "0_")

    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass
    cfg = jaws_site.config.Configuration(config_path)
    # Make sure that the _vars attribute in the interpolation object is "NONE" because
    # there was no env_override value set
    # print("Verifying that the _vars instance variable is NONE")
    assert cfg.config._vars is None
    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass

    # Check for case of env prefix set, but no variables match based on
    # https://code.jgi.doe.gov/advanced-analysis/jaws/jaws-central/-/issues/7
    cfg = jaws_site.config.Configuration(config_path, "ENV__")
    # Clear any pre-existing env vars
    try:
        del os.environ["JAWS_DB_HOST"]
        del os.environ["JAWS_DB_PASSWORD"]
    except KeyError:
        pass

    check_section("DB", expected_db_sections, cfg)
    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass

    os.environ["ENV__host2"] = "db_over.foobar.com"
    os.environ["ENV__password2"] = "over_password123"
    # print("Recreating new JAWSConfig object with env_override set to ENV__")
    cfg = jaws_site.config.Configuration(config_path, "ENV__")
    # print("Verifying that _vars attribute has picked up environment variables")
    assert cfg.config._vars is not None
    assert len(cfg.config._vars) == 2
    assert cfg.config._vars["host2"] == "db_over.foobar.com"
    assert cfg.config._vars["password2"] == "over_password123"

    # print("Retesting with env_prefix set")
    check_section("DB", expected_db_sections_env_override, cfg)

    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass
    # print("Testing with env_prefix set and section names")
    os.environ["JAWS_DB_HOST"] = "db.foobar.com"
    os.environ["JAWS_DB_PASSWORD"] = "password"
    os.environ["ENV__DB_host"] = "db-host.foo.com"
    os.environ["ENV__DB_test"] = "testing"

    # Verify we get a KeyError when mixing section names and non-section
    with pytest.raises(KeyError):
        jaws_site.config.Configuration(config_path, "ENV__")
    del os.environ["ENV__host2"]
    del os.environ["ENV__password2"]

    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass
    # Check [DB]host setting after clearing errors
    cfg = jaws_site.config.Configuration(config_path, "ENV__")
    assert "DB" in cfg.config._vars
    assert cfg.config["DB"]["host"] == os.environ["ENV__DB_host"]
    assert cfg.config["DB"]["test"] == os.environ["ENV__DB_test"]
    assert cfg.get_section("DB")["host"] == cfg.config["DB"]["host"]

    check_section("DB", expected_db_sections_env_override2, cfg)

    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass

    # Test with get_section() with env_prefix set
    os.environ["ENV__DB_host"] = "mysql.db.host"
    os.environ["ENV__DB_password"] = "mysqlpassword"
    cfg = jaws_site.config.Configuration(config_path, "ENV__")

    db_conf = cfg.get_section("DB")
    assert db_conf["host"] == os.environ["ENV__DB_host"]
    assert db_conf["password"] == os.environ["ENV__DB_password"]



@pytest.mark.parametrize(
    "section, key",
    [
        ("RMQ", "user"),
    ],
)
def test_get(partial_config, section, key):
    # call destructor to remove old references
    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass

    config_path = partial_config
    cfg = jaws_site.config.Configuration(config_path)
    assert cfg.config._vars is None
    assert cfg.get(section, key) == "bugs_bunny"


@pytest.mark.parametrize(
    "section, key",
    [
        ("RMQ_xxx", "user_xxx"),
    ],
)
def test_get_wrong_key(partial_config, section, key):
    # call destructor to remove old references
    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass

    config_path = partial_config
    cfg = jaws_site.config.Configuration(config_path)
    assert cfg.config._vars is None
    with pytest.raises(ConfigurationError):
        cfg.get(section, key)


@pytest.mark.parametrize(
    "section, key",
    [
        ("RMQ", "user"),
    ],
)
def test_get_exception(config_file_wrong, section, key):
    # call destructor to remove old references
    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass

    config_path = config_file_wrong
    with pytest.raises(ValueError):
        jaws_site.config.Configuration(config_path)


@pytest.mark.parametrize(
    "section",
    [
        ("RMQ"),
    ],
)
def test_get_section(partial_config, section):
    # call destructor to remove old references
    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass

    config_path = partial_config
    cfg = jaws_site.config.Configuration(config_path)
    assert cfg.config._vars is None
    # print(cfg.get_section(section))
    assert cfg.get_section(section) == {
        "host": "https://rmq.nersc.gov",
        "port": "5672",
        "user": "bugs_bunny",
        "password": "xqweasdasa",
        "num_threads": "5",
        "max_retries": "3",
        "vhost": "jaws_test",
    }


def test_file_not_found_config(file_not_found_config):
    try:
        jaws_site.config.Configuration._destructor()

    except Exception:
        pass
    config_path = file_not_found_config
    with pytest.raises(FileNotFoundError):
        jaws_site.config.Configuration(config_path)


def test_get_site_config(config_file):
    # Because Configuration is a singleton, we call a destructor method to
    # remove any old reference.
    try:
        jaws_site.config.Configuration._destructor()
    except Exception:
        pass
    config_path = config_file
    cfg = jaws_site.config.Configuration(config_path)

    result = cfg.get_site_config()
    for key in [
        "max_ram_gb",
        "inputs_dir",
        "access_group",
        "globus_host_path",
        "globus_endpoint",
    ]:
        assert key in result
