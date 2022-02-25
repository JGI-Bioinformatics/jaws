import logging
from datetime import datetime
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
        else:
            if not self.run.model:
                msg = f"Run info not found for run_id={run_id}"
                logger.error(msg)
                raise RunNotFoundError(msg)
            else:
                self.task = tasks.TaskLog(self.session, run_id=run_id)

    def get_result(self) -> None:
        """Get the result status for the given run_id provided during object instantiation."""

        run_logs = runs.get_run_status_logs(self.session, self.run_id)
        result = None
        for row in run_logs:
            if 'status_to' in row and row['status_to'] in ('succeeded', 'failed'):
                result = row['status_to']
                break
        return result

    def run_info(self) -> None:
        """
        Given a SQLAlchemy model for a Run, create a dict with the desired fields.
        :param run: Run object
        :type run: model
        :param is_admin: True if current user is an administrator
        :type is_admin: bool
        :param verbose: True if all fields desired
        :type verbose: bool
        :return: selected fields
        :rtype: dict
        """

        info = {
            "run_id": self.run.model.id,
            "user_id": self.run.model.user_id,
            "email": self.run.model.email,
            "submitted": self.run.model.submitted.strftime("%Y-%m-%d %H:%M:%S"),
            "updated": self.run.model.updated.strftime("%Y-%m-%d %H:%M:%S"),
            "status": self.run.model.status,
            "status_detail": jaws_constants.task_status_msg.get(self.run.model.status, ""),
            "result": self.get_result(),
            "site_id": config.conf.get("SITE", "id"),
        }
        return info

    def task_summary(self) -> None:
        """Get task summary info for the run based on the run_id provided during object instantiation."""

        try:
            entries = self.task.task_summary()
        except Exception as err:
            logger.error(f"Failed to get task_summary for run_id={self.run_id}: {err}")
            return {}

        infos = {}
        for entry in entries:
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

        return cromwell_metadata

    def create_doc(self) -> None:
        """Create a json document of the run info for the run_id provided during object instantiation."""

        run_status = self.run_info()
        if not run_status:
            return {}

        task_summary = self.task_summary()
        task_status = self.task_status()
        cromwell_metadata = self.cromwell_metadata()

        doc = {
            'run_id': run_status['run_id'],
            'workflow_name': cromwell_metadata.get('workflowName'),
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

    @staticmethod
    def _convert_date_to_seconds(time_string, format="%H:%M:%S") -> int:
        """Convert dates in string format to datetime object."""
        try:
            date_time = datetime.strptime(time_string, "%H:%M:%S")
        except ValueError:
            seconds = 0
        else:
            time_delta = date_time - datetime(1900, 1, 1)
            seconds = time_delta.total_seconds()
        return seconds

    @staticmethod
    def _mutate_run_dates(jsondata: dict) -> None:
        """Convert dates in string format to datetime object."""
        jsondata['submitted'] = datetime.fromisoformat(jsondata['submitted'])
        jsondata['updated'] = datetime.fromisoformat(jsondata['updated'])
        jsondata['runtime_sec'] = (jsondata['updated'] - jsondata['submitted']).total_seconds()

    @staticmethod
    def _mutate_task_status_dates(jsondata: dict) -> None:
        """Convert dates in string format to datetime object."""
        for task in jsondata:
            entries = jsondata[task]
            if 'timestamp' in entries and entries['timestamp']:
                entries['timestamp'] = datetime.fromisoformat(entries['timestamp'])

    def _mutate_task_summary_dates(self, jsondata: dict) -> None:
        """Convert dates in string format to datetime object."""
        for task in jsondata:
            entries = jsondata[task]
            if 'queued' in entries and entries['queued']:
                entries['queued'] = datetime.fromisoformat(entries['queued'])
            if 'queue_wait' in entries and entries['queue_wait']:
                entries['queue_wait_sec'] = self._convert_date_to_seconds(entries['queue_wait'])
            if 'runtime' in entries and entries['runtime']:
                entries['runtime_sec'] = self._convert_date_to_seconds(entries['runtime'])
            if 'maxtime' in entries and entries['maxtime']:
                entries['maxtime_sec'] = self._convert_date_to_seconds(entries['maxtime'])


class RPC_ES:
    """Class to create a RMQ connection for sending json payload to a queue where logstash will retrieve
    the payload and update elasticsearch with the doc.
    """
    def __init__(self, logger: logging) -> None:
        """Creates a RMQ connection based on the connection parameters specified in the config entry
        DASHBOARD_RPC_CLIENT. The config entries must have the following keys in the following example:

        [DASHBOARD_RPC_CLIENT]
        user: RabbitMQ user
        password: RabbitMQ password
        host: RabbitMQ host
        port: RabbitMQ port
        vhost: vhost name
        queue: queue name
        """

        # Get config paramters for rabbitmq connection in the config file.
        rpc_params = config.conf.get_section("DASHBOARD_RPC_CLIENT")

        self.rpc = rpc_es.RPCRequest(rpc_params, logger)

    def send_rpc_run_info(self, payload: dict) -> Tuple[dict, int]:
        """Sends request to RabbitMQ/RPC and wait for response. If response fails, return non-zero status_code.

        :param entries: A dictionary containing the rabbitmq connection information.
        :type entries: dict
        :param try_count: for logging number of tries only.
        :type try_count: int
        :return jsondata: dictionary of the results returned from RPC request.
        :rtype jsondata: dictionary
        :return status_code: zero if successful, non-zero if connection fails or return json contains error msg.
        :rtype status_code: int
        """

        try:
            jsondata = self.rpc.request(payload)
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
