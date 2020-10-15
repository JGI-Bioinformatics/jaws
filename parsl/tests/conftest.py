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

[SITE_RPC_CLIENT]
user = jaws
password = p4s5w0rd
host = rpc.server.com
vhost = j4w5
port = 56789
queue = high-prio
"""
    cfg.write_text(content)
    return cfg.as_posix()

