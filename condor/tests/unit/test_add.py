import pytest
import jaws_condor.add
from unittest import mock


def test_start_file_logger(tmp_path):
    tmp_log = tmp_path / "test.txt"
    jaws_condor.add.start_file_logger(tmp_log)
