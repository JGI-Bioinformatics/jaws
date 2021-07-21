import json
from datetime import datetime, timezone
import pytz
import os
import copy


class Runs:
    """
    A set of Runs.
    """

    def __init__(self, runs, params={}):
        self.runs = runs
        self.timezone = params.get("timezone", None)
        if self.timezone is None:
            self.timezone = os.environ.get("JAWS_TZ", None)

    def _utc_to_local(self, utc_datetime):
        """Convert UTC time to the local time zone. This should handle daylight savings.
           Param:: utc_datetime: a string of date and time "2021-07-06 11:15:17".
        """
        # The timezone can be overwritten with a environmental variable; this is useful if the
        # system time does not reflect the local time zone.
        # JAWS_TZ should be set to a timezone in a similar format to 'US/Pacific'
        local_tz = os.environ.get("JAWS_TZ", None)
        local_tz_obj = ''
        if local_tz is None:
            local_tz_obj = datetime.now().astimezone().tzinfo
        else:
            local_tz_obj = pytz.timezone(local_tz)

        fmt = "%Y-%m-%d %H:%M:%S"
        datetime_obj = datetime.strptime(utc_datetime, fmt)
        return datetime_obj.replace(tzinfo=timezone.utc).astimezone(tz=local_tz_obj).strftime(fmt)
