"""
This module provides a collection of functions called by rpc_operations.py.
"""

import json
import logging
import os
import re
import subprocess
from datetime import datetime, timedelta


def _sbatch_test_only():
    """
    Run shim for "sbatch --test-only" to get estimated queue-wait times from scheduler.
    The shim may not exist, depending on the backend being used (e.g. does not exist if using AWS batch).
    The outputs from the slurm command should look something like a json
    {
      "small":"",
      "medium":"sbatch: Job 1239738 to start at 2023-03-29T13:28:43 (...) nid00311 in partition genepool_shared",
      "large":"",
      "xlarge":"sbatch: Job 3909868 to start at 2023-03-30T04:20:01 (...) exvivo014 in partition exvivo"
    }
    """
    JAWS_BIN_DIR = os.getenv("JAWS_BIN_DIR")
    queue_wait_shim = os.path.join(JAWS_BIN_DIR, "queue-wait.sh")
    if not os.path.exists(queue_wait_shim):
        raise Exception(f"The path to the shim doesn't exist: {queue_wait_shim}")
    results = subprocess.run(["bash", queue_wait_shim], capture_output=True, text=True)
    return json.loads(results.stdout)


def check_queue_wait(logger):
    """
    Check the queue wait times for all the sites by using the "sbatch --test-only" command which returns an estimated
    start time for a job. The commands used here are asking for the same resources as
    condor asks for when creating slurm pools.
    Each site has a combination of small, medium, large, xlarge pools.

    :param logger: The logger passed from rpc_operations.py
    :type params: object
    :return: estimated start times for sm, med, lg & xlg pools in the format: 2023-03-29T10:40:28.
    :rtype: json
    """
    try:
        dictData = _sbatch_test_only()
    except Exception as error:
        logger.error(error)
        return

    now = datetime.now()
    result = {}
    pattern = re.compile(r"to\sstart\sat\s(\S+)")
    for (key, value) in dictData.items():
        if value:
            match = re.search(pattern, value)
            if match:
                hit = match.group(1)
                estimate = datetime.strptime(hit, "%Y-%m-%dT%H:%M:%S")
                sec = int(timedelta.total_seconds(estimate - now))
                result[key] = sec
    return result


if __name__ == "__main__":
    logger = logging.getLogger(__package__)
    check_queue_wait(logger)
