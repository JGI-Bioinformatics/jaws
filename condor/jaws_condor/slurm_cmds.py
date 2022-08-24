from jaws_condor.cmd_utils import run_sh_command
import logging
import pandas as pd


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
    def columns(self) -> list:
        """
        Returns the cloumns gotten from the 'squeue' command
        """
        return self.squeue_columns

    def call_slurm_command(self, command):
        """
        Genereic function to call a command and raise an error
        """
        stdout, stderr, returncode = run_sh_command(command, log=logger, show_stdout=False)

        # If we have an error return a dictionary with 0 for each type and state
        if returncode != 0:
            logger.error(f"Bad Return code from: {command}")
            raise SlurmCmdFailed(f"'{command}' command failed with {stderr}")

        return stdout, stderr, returncode

    def squeue(self) -> dict:
        """
        Runs the squeue command and returns a dictionary
        compatible with a pandas dataframe
        """
        command = f"{self.squeue_cmd}  -u {self.user_name} {str(self.extra_args)}"
        try:
            stdout, stderr, returncode = self.call_slurm_command(command)
        except SlurmCmdFailed:
            logging.error('squeue failed')
            stdout = ""

        # Create an empty dataframe with the right columns if nothing was returned or an error
        if stdout == "":
            return pd.DataFrame([], columns=[*self.columns, 'TIME_SEC']).to_dict()

        # Gets jobs from output by splitting on new lines
        jobs = [job.split() for job in stdout.split("\n")]
        try:
            df = pd.DataFrame(jobs, columns=self.columns)
        except AssertionError:
            logging.error(f'{jobs} mismatch {self.columns}')
            return pd.DataFrame([], columns=[*self.columns, 'TIME_SEC']).to_dict()
        # Drops rows if they have nan values
        df = df.dropna(axis=0)

        # Adds a new column of time in seconds left to run
        df['TIME_SEC'] = df["TIME_LEFT"].apply(slurm_time_to_sec)

        return df.to_dict()

    def scancel(self, job_id: int = 0, cluster: str = ""):

        cluster = "" if cluster is None else cluster
        if cluster != "":
            cluster = f"-M {cluster}"

        command = f"scancel {cluster} {job_id}"
        try:
            stdout, stderr, returncode = self.call_slurm_command(command)
            return {'stdout': stdout, 'stderr': stderr, 'returncode': returncode}
        except SlurmCmdFailed:
            logging.error('scancel failed')
            return {'stdout': "failed", 'stderr': "failed", 'returncode': 1}

    def sbatch(self, compute_type: str = "medium", cluster: str = ""):

        cluster = "" if cluster is None else cluster
        if cluster != "":
            cluster = f"-M {cluster}"

        command = f"sbatch --parsable {cluster} {self.script_path}/condor_worker_{compute_type}.job"
        try:
            stdout, stderr, returncode = self.call_slurm_command(command)
            return {'stdout': stdout, 'stderr': stderr, 'returncode': returncode}
        except SlurmCmdFailed:
            logging.error('sbatch failed')
            return {'stdout': "failed", 'stderr': "failed", 'returncode': 1}


def slurm_time_to_sec(time_str):
    # split off days first
    time_str = time_str.split("-")
    if len(time_str) > 1:
        days = int(time_str[0])
    else:
        days = 0

    # split time_str into HH:MM:SS
    time_str = time_str[-1].split(":")

    time_str_bits = {0: 1, 1: 60, 2: 60*60, 3: 60*60*24}
    total = 0
    # Run in reverse becasue we will always
    # have seconds and not always hours
    # 0 -> sec
    # 1 -> min
    # 2 -> hrs.
    # 3 -> days.
    for i, t in enumerate(time_str[::-1]):
        total += (time_str_bits[i]*int(t))

    if days > 0:
        total += (time_str_bits[3]*int(days))

    # Return total seconds
    return total
