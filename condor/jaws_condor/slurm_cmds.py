# from jaws_condor import config
from jaws_condor.cmd_utils import run_sh_command
import logging
import pandas as pd
import time

logger = logging.getLogger(__package__)


class SlurmError(Exception):
    pass


class SlurmCmdFailed(SlurmError):
    pass


class Slurm:
    def __init__(self, **kwargs):
        pass

    def squeue(self):
        pass

    def scancel(self):
        pass

    def sbatch(self):
        pass
