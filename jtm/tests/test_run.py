# test_run.py
from jaws_jtm.lib.run import run, make_dir, rm_dir


def test_run():
    assert run(["rm", "-rf", "./unittest"], dry_run=False) == 0


def test_make_dir():
    assert make_dir("./unittest", dry_run=False) == 0


def test_rm_dir():
    assert rm_dir("./unittest", dry_run=False) == 0
