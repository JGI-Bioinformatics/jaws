"""
This module provides a collection of functions called by rpc_operations.py.
"""


def status(logger):
    """
    Check the status of all the sites by using the "sbatch --test-only" command which returns an estimatated
    start time for a job. The jobs are asking for the same resources as condor asks for when creating slurm pools.
    Each site has a combination of  small, medium, large, xlarge pools.
    """

    logger.debug("Getting status of ")
    # test integration shim. Not full path.  Use jaws_bin_dir.

    path = "$JAWS_BIN_DIR/slurm-scripts.sh"

    return 
