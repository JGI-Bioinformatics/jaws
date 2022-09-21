from jaws_site import log
import pytest


def test_setup_logger(tmp_path):
    name = "testlog"
    file_name = tmp_path / "test.log"
    ret = log.setup_logger(name, file_name)
    assert ret is not None

    # Test invalid log level
    with pytest.raises(ValueError):
        log.setup_logger(name, file_name, log_level="test")

    # Test non-log file
    name = "testlog"
    file_name = None
    ret = log.setup_logger(name, file_name)
    assert "test.log" in ret.handlers[1].baseFilename
