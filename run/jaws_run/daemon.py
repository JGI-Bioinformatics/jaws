"""
JAWS Daemon process periodically checks on runs and performs actions to usher
them to the next state.
"""

import logging
import schedule
import time
from sqlalchemy.exc import SQLAlchemyError
from jaws_run import run, db


logger = logging.getLogger(__package__)


class DatabaseError(Exception):
    pass


active_run_states = [
    "uploading",
    "upload complete",
    "submitted",
    "queued",
    "running",
    "succeeded",
    "failed",
    "downloading",
]


class Daemon:
    """
    The daemon periodically checks on active runs, which may usher them to the next state.
    """

    def __init__(self):
        """
        Init daemon with schedule
        """
        schedule.every(10).seconds.do(self.check_active_runs)

    def start_daemon(self):
        """
        Start the scheduled loop.
        """
        while True:
            schedule.run_pending()
            time.sleep(1)

    def check_active_runs(self):
        """
        Check on current status of active runs.  This is run as a new thread.
        """
        # since this is a new thread, must get new session
        session = db.Session()

        # get list of active runs
        try:
            active_runs = (
                db.session.query(db.Run)
                .filter(db.Run.status.in_(active_run_states))
                .all()
            )
        except SQLAlchemyError as error:
            logger.exception(f"Error selecting active runs: {error}")
            raise DatabaseError(f"{error}")

        # init Run objects and have them check and update their status
        for row in active_runs:
            active_run = run.Run(row.id, session)
            active_run.check_status()

        session.close()
