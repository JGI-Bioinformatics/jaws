# test_run.py
import sys
import subprocess
from jaws_jtm.lib.run import run, make_dir, rm_dir, back_ticks


def test_run():
    assert run(["rm", "-rf", "./unittest"], dryRun=False) == 0


def test_make_dir():
    assert make_dir("./unittest", dryRun=False) == 0


def test_rm_dir():
    assert rm_dir("./unittest", dryRun=False) == 0


def test_back_ticks():
    cmd = "free"
    try:
        std_out = back_ticks(cmd, shell=True)
    except subprocess.CalledProcessError as msg:
        sys.stderr.write("Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode))
    assert std_out is not None
