import json
import os
import copy
from jaws_client.model import Model


class Runs(Model):
    """
    A set of Runs.
    """

    def __init__(self, runs, local_tz: None:
        self.runs = runs
        self.local_tz = local_tz
        if self.timezone is None:
            self.timezone = os.environ.get("JAWS_TZ", None)

    def output(self):
        runs = copy.deepcopy(self.runs)
        for run in runs:
            run[1] = _self.utc_to_local(run[1])
        return runs
