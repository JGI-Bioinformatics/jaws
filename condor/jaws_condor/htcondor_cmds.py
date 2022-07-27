# from jaws_condor import config
from jaws_condor.cmd_utils import run_sh_command
import logging
import pandas as pd
import time

logger = logging.getLogger(__package__)


class HTCondorError(Exception):
    pass


class CondorCmdFailed(HTCondorError):
    pass


class HTCondor:
    def __init__(self, **kwargs):
        # self.config = config.conf.get_section("POOL_MANAGER")  # TODO Get right section
        # self.site = self.config["SITE_ID"]
        self.columns = "ClusterId RequestMemory RequestCpus \
            CumulativeRemoteSysCpu CumulativeRemoteUserCpu \
            JobStatus NumShadowStarts JobRunCount RemoteHost JobStartDate QDate"

        # columns = self.config["CONDOR"]["condor_columns"]
        self.condor_q_cmd = f"condor_q -allusers -af {self.columns}"

        self.get_idle = 'condor_status -const "TotalSlots == 1" -af Machine'

    def condor_q(self):
        # Call the command
        condor_jobs = self.call_condor_command(self.condor_q_cmd)
        # get serialized and processed dataframe
        dataframe = process_condor_q(condor_jobs, self.condor_columns())
        return dataframe

    def condor_idle(self):
        idle_list = self.call_condor_command(self.get_idle)
        idle = [item for sublist in idle_list for item in sublist]
        return idle

    def condor_columns(self):
        # Splits the string columns to be used as header of dataframe
        return self.columns.split()

    def call_condor_command(self, condor_cmd):
        # Call the condor_q with the desired options
        stdout, stderr, returncode = run_sh_command(condor_cmd, log=logger, show_stdout=False)
        # The usual failures
        if returncode != 0:
            logger.critical(f"ERROR: failed to execute condor_q command: {condor_cmd}")
            raise CondorCmdFailed(f"condor '{condor_cmd}' command failed with {stderr}")

        return parse_condor_outputs(stdout)


def parse_condor_outputs(stdout):
    # split outputs by rows
    outputs = stdout.split("\n")
    # Split each row into columns
    condor_jobs = [job.split() for job in outputs]
    # Removes columns with no values (Usually the last column)
    condor_jobs = [q for q in condor_jobs if len(q) != 0]
    return condor_jobs


def process_condor_q(condor_jobs, columns):
    # if there are no jobs send back empty array of jobs
    if len(condor_jobs) == 0:
        return [{k: 0 for k in columns}]
    # Create a dataframe from the split outputs
    df = pd.DataFrame(condor_jobs, columns=columns)
    # Change the types
    df["JobStatus"] = df["JobStatus"].astype(int)
    df["RequestMemory"] = df["RequestMemory"].astype(float) / 1024
    df["RequestCpus"] = df["RequestCpus"].astype(float)
    df["CumulativeRemoteSysCpu"] = df["CumulativeRemoteSysCpu"].astype(float)
    df["CumulativeRemoteUserCpu"] = df["CumulativeRemoteUserCpu"].astype(float)

    now = int(time.time())
    df["JobStartDate"] = df["JobStartDate"].str.replace('undefined', str(now))
    df["total_running_time"] = now - df["JobStartDate"].astype(int)
    df["cpu_percentage"] = (((df['CumulativeRemoteSysCpu'] + df['CumulativeRemoteUserCpu']
                              ) / df['RequestCpus']) / df['total_running_time']) * 100
    df["total_q_time"] = df["JobStartDate"].astype(int) - df["QDate"].astype(int)

    # Returns a serilized/json version better for passing between objects
    # And for testing
    return df.to_dict('records')
