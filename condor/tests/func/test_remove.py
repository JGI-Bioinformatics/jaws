from deepdiff import DeepDiff
import os
from jaws_condor import remove
from jaws_condor.pool_manager_pandas import PoolManagerPandas


def test_number_of_workers_rm():
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

    old_workers = [
        0,  # 1
        0,  # 2
        0,  # 3
        0,  # 4
        -26,  # 5 -> 1 worker total, 30 running, go to minimum pool (4 nodes)
        -2,  # 6 -> 8 workers total, 10 running, remove 2
    ]

    poolman = PoolManagerPandas()
    for condor, slurm, workers in zip(condor_job_queue, slurm_workers, old_workers):
        _workers = poolman.need_cleanup(condor_job_queue=condor, slurm_workers=slurm, machine_size="medium")
        # print(_workers)
        assert workers == _workers


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
