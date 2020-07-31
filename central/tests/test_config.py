import jaws_central.config


def check_section(section_to_test, expected_entries, actual_config):
    for section, expected in expected_entries:
        assert actual_config.get(section_to_test, section) == expected


def test_check_all_values(config_file):

    config_path = config_file
    cfg = jaws_central.config.Configuration(config_path)

    expected_db_sections = [
        ("dialect", "mysql+mysqlconnector"),
        ("host", "db.foo.com"),
        ("port", "3306"),
        ("user", "jaws"),
        ("password", "passw0rd1"),
        ("db", "jaws"),
    ]

    expected_globus_sections = [
        ("client_id", "ZZZZ"),
    ]

    expected_site_lbnl_sections = [
        ("amqp_host", "rmq.jaws.gov"),
        ("amqp_user", "jaws"),
        ("amqp_password", "passw0rd2"),
        ("amqp_vhost", "jaws"),
        ("amqp_queue", "lbnl_rpc"),
        ("globus_endpoint", "XXXX"),
        ("globus_basepath", "\"/global/scratch/jaws\""),
        ("staging_subdir", "\"staging\""),
        ("max_ram_gb", "1024"),
    ]

    expected_site_nersc_sections = [
        ("amqp_host", "rmq.jaws.gov"),
        ("amqp_user", "jaws"),
        ("amqp_password", "passw0rd2"),
        ("amqp_vhost", "jaws"),
        ("amqp_queue", "nersc_rpc"),
        ("globus_endpoint", "YYYY"),
        ("globus_basepath", "\"/\""),
        ("staging_subdir", "\"/global/scratch/jaws/staging\""),
        ("max_ram_gb", "2048"),
    ]

    check_section("DB", expected_db_sections, cfg)
    check_section("GLOBUS", expected_globus_sections, cfg)
    check_section("SITE:LBNL", expected_site_lbnl_sections, cfg)
    check_section("SITE:NERSC", expected_site_nersc_sections, cfg)


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
        ("client_id", "ZZZZ"),
    ]

    check_section("DB", expected_db_sections, cfg)
    check_section("GLOBUS", expected_globus_sections, cfg)
