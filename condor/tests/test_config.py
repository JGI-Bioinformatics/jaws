import jaws_condor.config


def test_get_site_id(config_file):
    config_path = config_file
    cfg = jaws_condor.config.Configuration(config_path)
    expected_site_id = "LOCAL"
    actual_site_id = cfg.get_site_id()
    assert actual_site_id == expected_site_id
