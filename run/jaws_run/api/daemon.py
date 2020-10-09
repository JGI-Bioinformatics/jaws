"""
JAWS Analysis Service API

This service stores persistent Run information in a db and interacts with Cromwell.
"""

import logging
import schedule
import time
from sqlalchemy.exc import SQLAlchemyError
from jaws_run import config, db
from jaws_run.cromwell import Cromwell
from jaws_run.api.run import Run


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


class DatabaseError(Exception):
    pass


class Daemon:
    """
    The daemon periodically checks on active runs, which may usher them to the next state.
    This is the only way Runs may change state.
    """

    # only Runs in these states have an arc to another state
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

    def __init__(self):
        """
        Init daemon with schedule
        """
        schedule.every(10).seconds.do(self.__check_active_runs)

    def start(self):
        """
        Start the scheduled loop.
        """
        while True:
            schedule.run_pending()
            time.sleep(1)

    def __check_active_runs(self):
        """
        Check on current status of active runs.  This is run as a new thread.
        """
        # since this is a new thread, must get new session
        session = db.Session()

        # get list of active runs
        try:
            active_runs = (
                db.session.query(db.Run)
                .filter(db.Run.status.in_(self.active_run_states))
                .all()
            )
        except SQLAlchemyError as error:
            logger.exception(f"Error selecting active runs: {error}")
            raise DatabaseError(f"{error}")

        # have each Run check/update their status
        for row in active_runs:
            run = Run(row.id, session)
            run.check_status()

        session.close()
