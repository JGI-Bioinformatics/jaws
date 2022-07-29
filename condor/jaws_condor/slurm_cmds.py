# from jaws_condor import config
from jaws_condor.cmd_utils import run_sh_command
import logging
# import pandas as pd
# import time

logger = logging.getLogger(__package__)


class SlurmError(Exception):
    pass


class SlurmCmdFailed(SlurmError):
    pass


class Slurm:
    def __init__(self, user_name: str = "jaws_jtm", extra_args: str = None, script_path: str = "", **kwargs):
        self.squeue_columns = ["JOBID", "PARTITION", "NAME", "STATE", "TIME_LEFT", "NODELIST"]
        self.squeue_cmd = 'squeue --format="%.12i %.10P %.50j %.10T %.10L %.12R" --noheader'
        self.user_name = user_name
        self.extra_args = extra_args
        self.script_path = script_path

    @property
    def columns(self):
        return self.squeue_columns

    def squeue(self):
        command = f"{self.squeue_cmd}  -u {self.user_name} {str(self.extra_args)}"
        stdout, stderr, returncode = run_sh_command(command, log=logger, show_stdout=False)

        # If we have an error return a dictionary with 0 for each type and state
        if returncode != 0:
            logger.error(f"No output from squeue: {command}")
            raise SlurmCmdFailed(f"'{command}' command failed with {stderr}")

        # Gets jobs from output by splitting on new lines
        jobs = [job.split() for job in stdout.split("\n")]

        return jobs

    def scancel(self, job_id: int = 0, cluster: str = None):

        if cluster is not None:
            cluster = f"-M {cluster}"
        else:
            cluster = ""

        command = f"scancel {cluster} {job_id}"
        stdout, stderr, returncode = run_sh_command(command, log=logger, show_stdout=False)

        # If we have an error return a dictionary with 0 for each type and state
        if returncode != 0:
            logger.error(f"No output from squeue: {command}")
            raise SlurmCmdFailed(f"'{command}' command failed with {stdout} {stderr}")

        return True

    def sbatch(self, compute_type: str = "medium", cluster: str = None):

        if cluster is not None:
            cluster = f"-M {cluster}"
        else:
            cluster = ""

        command = f"sbatch --parsable {cluster} {self.script_path}/condor_worker_{compute_type}.job"
        stdout, stderr, returncode = run_sh_command(command, log=logger, show_stdout=False)

        # If we have an error return a dictionary with 0 for each type and state
        if returncode != 0:
            logger.error(f"No output from squeue: {command}")
            raise SlurmCmdFailed(f"'{command}' command failed with {stderr}")

        return stdout
