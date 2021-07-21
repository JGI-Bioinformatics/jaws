from datetime import datetime, timezone
import pytz


class Model:
    """Model baseclass"""

    def _utc_to_local(self, utc_datetime, local_tz=None):
        """Convert UTC time to the local time zone. This should handle daylight savings.
           Param:: utc_datetime: a string of date and time "2021-07-06 11:15:17".
        """
        local_tz_obj = pytz.timezone(local_tz) if local_tz else datetime.now().astimezone().tzinfo
        fmt = "%Y-%m-%d %H:%M:%S"
        datetime_obj = datetime.strptime(utc_datetime, fmt)
        return datetime_obj.replace(tzinfo=timezone.utc).astimezone(tz=local_tz_obj).strftime(fmt)
