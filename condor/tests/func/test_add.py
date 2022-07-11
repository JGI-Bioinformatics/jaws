from deepdiff import DeepDiff
import os
import functools
from jaws_condor import add


def test_start_file_logger():
    test_filename = "./test_add_log.log"
    add.start_file_logger(test_filename)
    logger = add.logger
    test_log_message = "EXAMPLE WARNING LOG MESSAGE"
    logger.warning(test_log_message)
    logger.disabled = True
    with open(test_filename, "r") as fh:
        complete_log = fh.read()
    os.remove(test_filename)
    assert test_log_message in complete_log


def test_collect_condor_jobs():
    condor_q_out = """1      50 GB    1024.0 B  16
2      510.0 GB    1024.0 B  16
3      110.0 GB    1024.0 B  16
4      210.0 GB    1024.0 B  32
5      220.0 GB    1024.0 B  32
"""
    ram_range_list = ["0-100", "100-200", "200-500", "500-1000"]
    res_dict = add.collect_condor_jobs(condor_q_out, ram_range_list)
    res_expected = [
        [[1, 50.0, 1024.0, 16]],
        [[3, 110.0, 1024.0, 16]],
        [[4, 210.0, 1024.0, 32], [5, 220.0, 1024.0, 32]],
        [[2, 510.0, 1024.0, 16]]
    ]
    assert functools.reduce(lambda x, y: x and y, map(lambda p, q: p == q, res_dict, res_expected), True)


def test_calculate_node_needed():
    """
    $ condor_q -pr ../fmt_nobatch_id.cpf
    190     1024.0 KB    1024.0 B  1
    191     1024.0 KB    1024.0 B  1
    192      390.6 GB    1024.0 B  4
    193       97.7 GB    1024.0 B  10
    194       97.7 GB    1024.0 B  10
    195       97.7 GB    1024.0 B  10
    196       97.7 GB    1024.0 B  10

    ==> idle_list
    [
    [],
    [[190, 0.0009765625, 1024.0, 1], [191, 0.0009765625, 1024.0, 1], [193, 97.7, 1024.0, 10], [194, 97.7, 1024.0, 10]], [194, 97.7, 1024.0, 10]],  # noqa
    [],
    [[192, 390.6, 1024.0, 4]]
    ]

    ==> Resources per type
    types = ["small", "medium", "large", "xlarge"]
    ; unit GB
    mem = ["0-0", "0-364", "0-0", "364-1480"]
    cpu = [36, 36, 36, 36]

    ==> total number of nodes needed to sbatch
    max(sum(mem_requested) / each_max_mem, sum(cpu_requested) / each_max_cpu) +  1

    ==> final sbatch numbers
    [0, 2, 0, 1]
    """
    job_list = [[], [[190, 0.0009765625, 1024.0, 1], [191, 0.0009765625, 1024.0, 1], [193, 97.7, 1024.0, 10], [194, 97.7, 1024.0, 10], [195, 97.7, 1024.0, 10], [196, 97.7, 1024.0, 10]], [], [[192, 390.6, 1024.0, 4]]]  # noqa
    r_type = 1  # MEDIUM
    mem_range = ["0-0", "0-364", "0-0", "364-1480"]
    cpu_range = [36, 36, 36, 36]
    summed = [sum(x) for x in zip(*job_list[r_type])]
    assert [1159, 390.801953125, 6144.0, 42] == summed

    ret = add.calculate_node_needed(job_list[r_type], mem_range[r_type], cpu_range[r_type])
    assert ret == 2

    r_type = 3  # XLARGE
    ret = add.calculate_node_needed(job_list[r_type], mem_range[r_type], cpu_range[r_type])
    assert ret == 1
