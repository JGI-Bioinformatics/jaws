from deepdiff import DeepDiff
import os
from jaws_condor import remove


def test_start_file_logger():
    test_filename = "./test_remove_log.log"
    remove.start_file_logger(test_filename)
    logger = remove.logger
    test_log_message = "EXAMPLE WARNING LOG MESSAGE"
    logger.warning(test_log_message)
    logger.disabled = True
    with open(test_filename, "r") as fh:
        complete_log = fh.read()
    os.remove(test_filename)
    assert test_log_message in complete_log


def test_collect_condor_running_jobs():
    test_condor_q_out = """
190        1.0 TB    1024.0 MB  1
191     1024.0 KB    1024.0 MB  16
192      390.6 GB    1024.0 MB  4
193       97.7 GB    1024.0 MB  10
"""
    test_ram_range = ["0-0", "0-118", "0-0", "118-1450"]
    actual = remove.collect_condor_running_jobs(test_condor_q_out, test_ram_range)
    expected = [
        [],
        [["191", 0.0009765625, 1024.0, 16], ["193", 97.7, 1024.0, 10]],
        [],
        [["190", 1024.0, 1024.0, 1], ["192", 390.6, 1024.0, 4]],
    ]
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False
