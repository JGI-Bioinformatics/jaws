import json
from datetime import datetime, timezone
import pytz
import os
import copy


class RunLog:
    """
    A Run log.
    """
    def __init__(self, log: list):
        self.log = log

    def output(self, fmt = "text"):
        """
        Return log in desired format.
        """
        header = ["#STATUS_FROM", "STATUS_TO", "TIMESTAMP", "REASON"]
        if fmt == "json":
            return self.log
        elif fmt == "tab":
            table = "\t".join(header)
            for log_entry in self.log:
                table += "\t".join(log_entry)
            return table
        else:
            """Get the max length of element in every col and add padding (2)"""
            log = copy.deepcopy(self.log)
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
