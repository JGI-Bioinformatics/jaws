from deepdiff import DeepDiff
from jaws_condor import add


def test_collect_condor_jobs():
    test_condor_q_out = """
190        1.0 TB    1024.0 MB  1
191     1024.0 KB    1024.0 MB  16
192      390.6 GB    1024.0 MB  4
193       97.7 GB    1024.0 MB  10
"""
    test_ram_range = ["0-0", "0-118", "0-0", "118-1450"]
    actual = add.collect_condor_jobs(test_condor_q_out, test_ram_range)
    expected = [
        [],
        [[191, 0.0009765625, 1024.0, 16], [193, 97.7, 1024.0, 10]],
        [],
        [[190, 1024.0, 1024.0, 1], [192, 390.6, 1024.0, 4]],
    ]
    assert bool(DeepDiff(actual, expected, ignore_order=True)) is False


# def test_calculate_node_needed():
#     test_idle_list = [
#         [],
#         [[191, 0.0009765625, 1024.0, 16], [193, 97.7, 1024.0, 10]],
#         [],
#         [[190, 1024.0, 1024.0, 1], [192, 390.6, 1024.0, 4]],
#     ]
#     test_ram_range = "118-1450"
#     test_ncpu = 36
#     actual = add.calculate_node_needed(test_idle_list, test_ram_range, test_ncpu)
#     print(actual)
#     expected = 0
#     assert actual == expected
