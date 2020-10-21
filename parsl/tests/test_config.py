import jaws_parsl.config


def test_check_rmq_values(config_file):

    config_path = config_file
    cfg = jaws_parsl.config.Configuration(config_path)

    expected_rmq_parameters = [
        ("user", "j4w5"),
        ("password", "p455w0rd"),
        ("host", "rmq.server.com"),
        ("vhost", "jaws"),
        ("port", 5678)
    ]

    actual_parameters = cfg.get_rmq_params()

    for key, expected in expected_rmq_parameters:
        assert actual_parameters[key] == expected
