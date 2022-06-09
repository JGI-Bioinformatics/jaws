import datetime


def today(format: str = "Y%-%m-%d") -> str:
    """Return the current date in string format"""
    return to_str_format(datetime.datetime.now())


def subtract_days(date: str, num_days: int, format: str = "%Y-%m-%d") -> str:
    """Subtract num_days to input date and return new date string"""
    date = to_datetime_format(date, format)
    new_date = date - datetime.timedelta(days=num_days)
    return to_str_format(new_date)


def add_days(date: str, num_days: int, format: str = "%Y-%m-%d") -> str:
    """Add num_days to input date and return new date string"""
    date = to_datetime_format(date, format)
    new_date = date + datetime.timedelta(days=num_days)
    return to_str_format(new_date)


def to_str_format(date: datetime, format: str = "%Y-%m-%d") -> str:
    """Convert datetime object to string"""
    return date.strftime(format)


def to_datetime_format(date: str, format: str = "%Y-%m-%d") -> datetime:
    """Convert string to datetime object"""
    return datetime.datetime.strptime(date, format)


def get_seconds_between_dates(start_date: str, end_date: str, date_format: str = "%Y-%m-%d") -> int:
    """Given two date strings, return number of secs between two dates"""
    start_date = strip_timezone_from_date(start_date)
    end_date = strip_timezone_from_date(end_date)

    start_dt = datetime.datetime.strptime(start_date, date_format)
    end_dt = datetime.datetime.strptime(end_date, date_format)

    return (end_dt - start_dt).total_seconds()


def strip_timezone_from_date(date: str) -> str:
    """If date string contains timezone like 2022-05-12T17:54:27.778Z, remove timezone from string"""
    if "." in date:
        idx = date.index(".")
        if idx >= 0:
            date = date[0:idx]
    return date
