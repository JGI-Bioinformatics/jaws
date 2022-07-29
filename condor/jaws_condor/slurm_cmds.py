# from jaws_condor import config
# from jaws_condor.cmd_utils import run_sh_command
import logging
# import pandas as pd
# import time

logger = logging.getLogger(__package__)


class SlurmError(Exception):
    pass


class SlurmCmdFailed(SlurmError):
    pass


class Slurm:
    def __init__(self, **kwargs):
        self.squeue_columns = ["JOBID", "PARTITION", "NAME", "STATE", "TIME_LEFT", "NODELIST"]
        self.squeue_cmd = 'squeue --format = "%.12i %.10P %.50j %.10T %.10L %.12R" --noheader'
        pass

    @property
    def columns(self):
        return self.squeue_columns

    def squeue(self, user_name, options):
        # extra = f"-u {user_name} {options}"
        pass

    def scancel(self):
        pass

    def sbatch(self):
        pass
