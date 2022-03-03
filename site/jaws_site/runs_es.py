import logging
from typing import Tuple
from jaws_rpc import responses
from jaws_site import (
    config,
    jaws_constants,
    rpc_es,
    database,
    runs,
    tasks
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

    def __init__(self, run_id: int) -> None:
        """Initialize database for retrieving run info.

        :param run_id: jaws run_id
        :type run_id: int
        """

        self.session = database.Session()
        self.run_id = run_id

        try:
            self.run = runs.Run(self.session, run_id=run_id)
        except (runs.RunNotFound, runs.RunDbError) as err:
            msg = f"Failed to get run info for run_id={run_id}: {err}"
            logger.error(msg)
            raise RunNotFoundError(msg)

        self.task = tasks.TaskLog(self.session, run_id=run_id)

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
        """Get cromwell metadata for the run based on the run_id provided during object instantiation."""

        try:
            cromwell_metadata = self.run.metadata()
        except Exception as err:
            logger.error(f"Failed to get task_status for run_id={self.run_id}: {err}")
            return {}

        return cromwell_metadata or {}

    def create_doc(self) -> None:
        """Create a json document of the run info for the run_id provided during object instantiation."""

        run_status = self.run.get_run_metadata()
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
