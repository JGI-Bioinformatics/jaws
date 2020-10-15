import jaws_auth.config


def check_section(section_to_test, expected_entries, actual_config):
    for key, expected in expected_entries:
        assert actual_config.get(section_to_test, key) == expected


def test_check_all_values(config_file):

    config_path = config_file
    cfg = jaws_auth.config.Configuration(config_path)

    expected_db_parameters = [
        ("dialect", "mysql+mysqlconnector"),
        ("host", "db.foo.com"),
        ("port", "3306"),
        ("user", "jaws"),
        ("password", "passw0rd1"),
        ("db", "jaws"),
    ]

    expected_auth_parameters = [
        ("port", "3000"),
    ]

    check_section("DB", expected_db_parameters, cfg)
    check_section("AUTH", expected_auth_parameters, cfg)


def test_config_overwrite_partial_values(partial_config):

    # call destructor to remove old references
    jaws_auth.config.Configuration._destructor()

    config_path = partial_config
    cfg = jaws_auth.config.Configuration(config_path)

    expected_db_parameters = [
        ("dialect", "mysql+mysqlconnector"),
        ("host", "db.foo.com"),
        ("port", "3305"),
        ("user", "j4ws"),
        ("password", "p455w0rd1"),
        ("db", "jaws"),
    ]
    expected_auth_parameters = [
        ("port", "3000")
    ]

    check_section("DB", expected_db_parameters, cfg)
    check_section("AUTH", expected_auth_parameters, cfg)
