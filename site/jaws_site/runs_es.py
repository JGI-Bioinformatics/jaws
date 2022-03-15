import logging
from typing import Tuple, Callable
from jaws_rpc import responses
from jaws_site import (
    rpc_es,
    tasks,
    models,
    cromwell,
)
from jaws_rpc.rpc_client import (
    InvalidJsonResponse,
    ConfigurationError,
    ConnectionError
)

logger = logging.getLogger(__package__)


class RunNotFoundError(Exception):
    pass


class RunES:
    """Class to retrieve jaws run information for a given run_id and create a json document to
    insert into elasticsearch."""

    def __init__(self, session: Callable, model: models.Run = None, run_id: int = None) -> None:
        """Initialize database for retrieving run info.

        :param session: database.Session() which is a sqlalchemy session object.
        :type session: sqlalchemy/database.session object
        :param run_id: jaws run_id
        :type run_id: int
        """

        self.session = session
        self.run_id = run_id
        self.model = None

        if model:
            self.model = model
            self.run_id = model.id
        elif run_id:
            """Select run record from rdb"""
            try:
                self.model = self.session.query(models.Run).get(run_id)
            except IntegrityError as error:
                logger.warn(f"Run {run_id}: IntegrityError for model thrown: {error}")
                raise RunNotFound(f"Run {run_id} not found")
            except SQLAlchemyError as error:
                err_msg = f"Unable to select run, {run_id}: {error}"
                logger.error(err_msg)
                raise RunDbError(err_msg)

        if not self.model:
            logger.warn(f"Run {run_id}: model not found")
            raise RunNotFoundError(f"Run {run_id} not found")

        self.task = tasks.TaskLog(self.session, run_id=self.run_id)

    def task_summary(self) -> None:
        """Get task summary info for the run based on the run_id provided during object instantiation."""

        try:
            entries = self.task.task_summary()
        except Exception as err:
            logger.error(f"Failed to get task_summary for run_id={self.run_id}: {err}")
            return {}

        infos = {}
        for entry in entries:
            if len(entry) < 8:
                continue
            task_name = entry[0]
            infos[task_name] = {
                'cromwell_id': entry[1],
                'cached': entry[2],
                'result': entry[3],
                'queued': entry[4],
                'queue_wait': entry[5],
                'runtime': entry[6],
                'maxtime': entry[7],
            }
        return infos

    def task_status(self) -> None:
        """Get task status info for the run based on the run_id provided during object instantiation."""

        try:
            entries = self.task.task_status()
        except Exception as err:
            logger.error(f"Failed to get task_status for run_id={self.run_id}: {err}")
            return {}

        infos = {}
        for entry in entries:
            if len(entry) < 6:
                continue
            task_name = entry[0]
            infos[task_name] = {
                'cromwell_id': entry[1],
                'status': entry[3],
                'timestamp': entry[4],
                'reason': entry[5],
            }
        return infos

    def cromwell_metadata(self) -> None:
        """
        Get metadata from Cromwell and return it, if available.
        If the run hasn't been submitted to Cromwell yet, the result shall be None.
        """
        if self.model.cromwell_run_id:
            return cromwell.get_metadata(self.model.cromwell_run_id).data
        else:
            return None

        return cromwell_metadata or {}

    def create_doc(self) -> None:
        """Create a json document of the run info for the run_id provided during object instantiation."""

        try:
            run = runs.Run(self.session, run_id=self.run_id)
        except (runs.RunNotFound, runs.RunDbError) as err:
            msg = f"Failed to get run info for run_id={run_id}: {err}"
            logger.error(msg)
            raise RunNotFoundError(msg)

        run_status = run.get_run_metadata()
        if not run_status:
            return {}

        task_summary = self.task_summary()
        task_status = self.task_status()
        cromwell_metadata = self.cromwell_metadata()

        doc = {
            'run_id': run_status['run_id'],
            'workflow_name': cromwell_metadata.get('workflowName', 'unknown'),
        }
        doc.update(run_status)
        doc['tasks'] = []

        for task_name in task_status:
            status_entries = {}
            status_entries['name'] = task_name
            status_entries.update(task_status[task_name])
            if task_name in task_summary:
                status_entries.update(task_summary[task_name])
            doc['tasks'].append(status_entries)

        return doc


def send_rpc_run_metadata(rpc_client: rpc_es.RPCRequest, payload: dict) -> Tuple[dict, int]:
    """Sends request to RabbitMQ/RPC and wait for response. If response fails, return non-zero status_code.

    :param rpc_client: rpc object connected to a RMQ queue for publishing message.
    :type rpc_client: rpc_es.RPCRequest object.
    :param payload: json document to publish to RMQ queue.
    :type payload: dict
    :return jsondata: response from RMQ request.
    :rtype jsondata: dictionary
    :return status_code: zero if successful, non-zero if connection fails or return json contains error msg.
    :rtype status_code: int
    """

    try:
        jsondata = rpc_client.request(payload)
    except InvalidJsonResponse as err:
        msg = f"RPC request returned an invalid response: {err}"
        logger.debug(msg)
        jsondata = responses.failure(err)
    except ConfigurationError as err:
        msg = f"RPC request returned an invalid configuration error: {err}"
        logger.debug(msg)
        jsondata = responses.failure(err)
    except ConnectionError as err:
        msg = f"RPC request returned an invalid connection error: {err}"
        logger.debug(msg)
        jsondata = responses.failure(err)

    status_code = 0

    if jsondata and 'error' in jsondata:
        status_code = jsondata['error'].get('code', 500)

    return jsondata, status_code
