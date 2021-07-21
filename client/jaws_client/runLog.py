import json
import os
import copy
from jaws_client.model import Model


class RunLog(Model):
    """
    A Run log.
    """
    def __init__(self, log: list, local_tz: string = None):
        self.log = log
        self.local_tz = local_tz

    def output(self, fmt = "text", local_tz: str = None):
        """
        Return log in desired format.
        """
        log = copy.deepcopy(self.log)
        for row in log:
            row[2] = self._utc_to_local(row[2])
        header = ["#STATUS_FROM", "STATUS_TO", "TIMESTAMP", "REASON"]
        if fmt == "json":
            return log
        elif fmt == "tab":
            table = "\t".join(header)
            for log_entry in log:
                table += "\t".join(log_entry)
            return table
        else:
            """Get the max length of element in every col and add padding (2)"""
            log.insert(0, header)
            col_widths = []
            for idx in range(len(header)):
                col_widths.append(max(len(log_entry[idx]) for log_entry in log) + 2)
            table = ""
            for log_entry in log:
                table += "".join(
                    cell.ljust(col_widths[col_idx])
                    for col_idx, cell in enumerate(log_entry)
                )
            return table
