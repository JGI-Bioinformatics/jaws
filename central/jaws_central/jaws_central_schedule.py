#!/usr/bin/env python3

"""
JAWS Central needs to perform regularly scheduled tasks.
Each task is it's own function.  We are using this instead of cron, which isn't available on all systems and we want all systems to be maintained as nearly identically as possible.
"""

import sys
import schedule
import time
import os


def refresh_globus_tokens()
    """
    Periodically refresh all users' Globus tokens so that they don't expire.
    """
    # NYI



## SCHEDULE
#schedule.every(10).seconds.do(hello)
#schedule.every(1).minutes.do(taskX)
#schedule.every(1).hour.do(something_else)
#schedule.every(1).day.at("23:00").do(some_daily_task)
#schedule.every().sunday.at("4:00").do(some_weekly_task)

schedule.every().sunday.at("4:00").do(refresh_globus_tokens)


while True:
    schedule.run_pending()
    time.sleep(60)
