import pytest


@pytest.fixture
def script(tmp_path):
    scr = tmp_path / "script.sh"
    content = """#!/bin/bash
echo 'hello'
"""
    scr.write_text(content)
    return scr.as_posix()


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "parsl.ini"
    content = """[RMQ]
user = j4w5
password = p455w0rd
host = rmq.server.com
vhost = jaws
port = 5678
"""
    cfg.write_text(content)
    return cfg.as_posix()

