import os
import functools
from jaws_condor import add
from jaws_condor.pool_manager_pandas import PoolManagerPandas
from jaws_condor.htcondor_cmds import HTCondor
from jaws_condor.slurm_cmds import Slurm


def test_number_of_workers_add():
    condor_job_queue = [
        {  # 1
            'hold_and_impossible': 0,
            'idle_medium': 5,
            'running_medium': 60,
            'medium_cpu_needed': 2064.0,
            'medium_mem_needed': 3575.0,
            'idle_xlarge': 0,
            'running_xlarge': 0,
            'xlarge_cpu_needed': 0,
            'xlarge_mem_needed': 0
        },
        {  # 2
            'hold_and_impossible': 0,
            'idle_medium': 5,
            'running_medium': 60,
            'medium_cpu_needed': 2064.0,
            'medium_mem_needed': 3575.0,
            'idle_xlarge': 0,
            'running_xlarge': 0,
            'xlarge_cpu_needed': 0,
            'xlarge_mem_needed': 0
        },
        {  # 3
            'hold_and_impossible': 0,
            'idle_medium': 5,
            'running_medium': 60,
            'medium_cpu_needed': 2064.0,
            'medium_mem_needed': 3575.0,
            'idle_xlarge': 0,
            'running_xlarge': 0,
            'xlarge_cpu_needed': 0,
            'xlarge_mem_needed': 0
        },
        {  # 4
            'hold_and_impossible': 0,
            'idle_medium': 1,
            'running_medium': 0,
            'medium_cpu_needed': 4.0,
            'medium_mem_needed': 8.0,
            'idle_xlarge': 0,
            'running_xlarge': 0,
            'xlarge_cpu_needed': 0,
            'xlarge_mem_needed': 0
        },
        {  # 5
            'hold_and_impossible': 0,
            'idle_medium': 1,
            'running_medium': 0,
            'medium_cpu_needed': 4.0,
            'medium_mem_needed': 8.0,
            'idle_xlarge': 0,
            'running_xlarge': 0,
            'xlarge_cpu_needed': 0,
            'xlarge_mem_needed': 0
        },
        {  # 6
            'hold_and_impossible': 0,
            'idle_medium': 0,
            'running_medium': 0,
            'medium_cpu_needed': 0,
            'medium_mem_needed': 60*8,  # Need 8
            'idle_xlarge': 0,
            'running_xlarge': 0,
            'xlarge_cpu_needed': 0,
            'xlarge_mem_needed': 0
        }
    ]

    slurm_workers = [
        {  # 1
            'medium_pending': 0,
            'medium_running': 30,  # Max pool
            'xlarge_pending': 0,
            'xlarge_running': 0
        },
        {  # 2
            'medium_pending': 0,
            'medium_running': 10,
            'xlarge_pending': 0,
            'xlarge_running': 0
        },
        {  # 3
            'medium_pending': 10,
            'medium_running': 10,
            'xlarge_pending': 0,
            'xlarge_running': 0
        },
        {  # 4
            'medium_pending': 0,
            'medium_running': 4,  # Min pool
            'xlarge_pending': 0,
            'xlarge_running': 0
        },
        {  # 5
            'medium_pending': 0,
            'medium_running': 30,  # Max pool
            'xlarge_pending': 0,
            'xlarge_running': 0
        },
        {  # 6
            'medium_pending': 0,
            'medium_running': 10,
            'xlarge_pending': 0,
            'xlarge_running': 0
        }
    ]

    new_workers = [
        0,  # 1
        20,  # 2
        10,  # 3
        0,  # 4
        0,  # 5
        0,  # 6
    ]
    configs = {
        "compute_types": [
            "medium",
            "xlarge"
        ],
        "user_name": "jaws_jtm",
        "squeue_args": "--clusters=all -p genepool,genepool_shared,exvivo,exvivo_shared",
        "min_pool": {
            "medium": 4,
            "xlarge": 0
        },
        "max_pool": {
            "medium": 30,
            "xlarge": 10
        },
        "worker_sizes": {
            "medium_cpu": 64,
            "medium_mem": 120,
            "xlarge_cpu": 72,
            "xlarge_mem": 1500
        },
        "cpu_bins": [
            0,
            64,
            72,
            10000
        ],
        "mem_bins": [
            0,
            120,
            1500,
            10000
        ],
        "labels": [
            "medium",
            "xlarge",
            "over"
        ]
    }
    wanted_columns = "ClusterId RequestMemory RequestCpus CumulativeRemoteSysCpu CumulativeRemoteUserCpu JobStatus NumShadowStarts JobRunCount RemoteHost JobStartDate QDate"  # noqa
    poolman = PoolManagerPandas(condor_provider=HTCondor(columns=wanted_columns),
                                slurm_provider=Slurm(), configs=configs)

    for condor, slurm, workers in zip(condor_job_queue, slurm_workers, new_workers):
        _workers = poolman.need_new_nodes(condor_job_queue=condor, slurm_workers=slurm, machine_size="medium")
        # print(_workers)
        assert workers == _workers


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
