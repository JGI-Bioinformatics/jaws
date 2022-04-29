import jaws_condor.add
import functools


def test_collect_condor_jobs():
    condor_q_out = """1      50 GB    1024.0 B  16
2      510.0 GB    1024.0 B  16
3      110.0 GB    1024.0 B  16
4      210.0 GB    1024.0 B  32
5      220.0 GB    1024.0 B  32
"""
    ram_range_list = ["0-100", "100-200", "200-500", "500-1000"]
    res_dict = jaws_condor.add.collect_condor_jobs(condor_q_out, ram_range_list)
    res_expected = [
        [[1, 50.0, 1024.0, 16]],
        [[3, 110.0, 1024.0, 16]],
        [[4, 210.0, 1024.0, 32], [5, 220.0, 1024.0, 32]],
        [[2, 510.0, 1024.0, 16]]
    ]
    assert functools.reduce(lambda x, y: x and y, map(lambda p, q: p == q, res_dict, res_expected), True)
