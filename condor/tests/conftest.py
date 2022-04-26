import pytest


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "jaws_condor.ini"
    content = """[SITE]
site_id = LOCAL

"""
    cfg.write_text(content)
    return cfg.as_posix()
