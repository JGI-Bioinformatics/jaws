from jaws_central.models_fsa import db, Run_Log

class RunLog:

    def __init__(self, session, run_id):
        self.session = session
        self.id = run_id

    def run_log(self):
        """
        Retrieve complete log of a Run's state transitions.

        :return: Table of log entries
        :rtype: list
        """
        logger.info(f"User {self.user.id}: Get log of Run {self.id}")
        try:
            query = (
                db.session.query(Run_Log)
                .filter_by(run_id=run_id)
                .order_by(Run_Log.timestamp)
            )
        except SQLAlchemyError as error:
            logger.exception(f"Error selecting from run_logs: {error}")
            abort(500, f"Db error; {error}")
        table = []
        for log in query:
            reason = log.reason if log.reason else ""
            row = [
                log.status_from,
                log.status_to,
                log.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                reason,
            ]
            table.append(row)
        return table, 200

    def add_log(self, status_from, status_to, timestamp, reason=None): 
        log = Run_Log(
            run_id=run.id,
            status_from=status_from,
            status_to=run.status,
            timestamp=run.updated,
            reason=reason,
        )
        try:
            db.session.add(log)
            db.session.commit()
        except Exception as error:
            db.session.rollback()
            logger.error(
                f"Error while adding run log entry to cancel run {run.id}: {error}"
            )
