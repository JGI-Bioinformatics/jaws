"""
This module provides a collection of functions called by rpc_operations.py.
"""

import os
import re
import subprocess
import json
import logging


def check_queue_wait(logger):
    """
    Check the queue wait times for all the sites by using the "sbatch --test-only" command which returns an estimatated
    start time for a job. The commands used here are asking for the same resources as
    condor asks for when creating slurm pools.
    Each site has a combination of small, medium, large, xlarge pools.

    :param logger: The logger passed from rpc_operations.py
    :type params: object
    :return: estimated start times for sm, med, lg & xlg pools in the format: 2023-03-29T10:40:28.
    :rtype: json
    """

    JAWS_BIN_DIR = os.getenv('JAWS_BIN_DIR')
    queue_wait_shim = os.path.join(JAWS_BIN_DIR, "queue-wait.sh")
    if not os.path.exists(queue_wait_shim):
        raise Exception(f"The path to the shim doesn't exist: {queue_wait_shim}")

    try:
        results = subprocess.run(["bash", queue_wait_shim], capture_output=True, text=True)
    except Exception as error:
        logger.error(error)
        return

    """
    The outputs from the slurm command should look something like a json
    {
      "small":"",
      "medium":"sbatch: Job 1239738 to start at 2023-03-29T13:28:43 (...) nid00311 in partition genepool_shared",
      "large":"",
      "xlarge":"sbatch: Job 3909868 to start at 2023-03-30T04:20:01 (...) exvivo014 in partition exvivo"
    }
    """

    dictData = json.loads(results.stdout)

    times = {}
    for (k, v) in dictData.items():
        if v:
            try:
                x = re.search(r"to\sstart\sat\s(\S+)", v)
                hit = x.group(1)
            except Exception as error:
                logger.error(f"failed to parse start time from the slurm command: {error}")

            times[k] = hit

    return times


if __name__ == "__main__":
    logger = logging.getLogger(__package__)
    check_queue_wait(logger)
