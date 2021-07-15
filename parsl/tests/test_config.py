import jaws_parsl.config


def test_check_mp_values(config_file):

    config_path = config_file
    cfg = jaws_parsl.config.Configuration(config_path)

    expected_mp_parameters = [
        ("password", "p455w0rd"),
        ("host", "mp.server.com"),
        ("port", 5678)
    ]

    actual_parameters = cfg.get_mp_params()

    for key, expected in expected_mp_parameters:
        assert actual_parameters[key] == expected


def test_check_rpc_values(config_file):

    config_path = config_file
    cfg = jaws_parsl.config.Configuration(config_path)

    expected_rpc_parameters = [
        ("user", "jaws"),
        ("password", "p4s5w0rd"),
        ("host", "rpc.server.com"),
        ("vhost", "j4w5"),
        ("port", 56789),
        ("queue", "high-prio")
    ]

    actual_parameters = cfg.get_rpc_params()

    for key, expected in expected_rpc_parameters:
        assert actual_parameters[key] == expected
