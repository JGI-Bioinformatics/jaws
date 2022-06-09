import logging
import elasticsearch
import elastic_transport
from typing import List
from jaws_central import config

logger = logging.getLogger(__package__)


class ElasticsearchError(elasticsearch.TransportError):
    pass


class ConfigurationError(Exception):
    pass


class ESClient:
    def __init__(self):
        es_params = config.conf.get_section("ELASTIC_SEARCH")
        es_indices = config.conf.get_section("RUN_METRICS_ELASTIC_INDEX")

        self._check_for_required_keys(es_params, ["host", "port", "api_key"])
        self._check_for_required_keys(es_indices, ["runs", "performance_metrics"])

        self.es_client = self.es_client(es_params["host"], es_params["port"], es_params["api_key"])
        self.es_index = es_indices

    @staticmethod
    def _check_for_required_keys(params: dict, required_params: list) -> None:
        """Given a dict and a list of keys to check within the dict, raise ConfigurationError if required
        key is missing in the dict.
        """
        if not all(key in params for key in required_params):
            raise ConfigurationError(f"{required_params} required")

    @staticmethod
    def es_client(host: str, port: int, api_key: str) -> elasticsearch.Elasticsearch:
        """Connect to elasticsearch and return elasticsearch client object.

        :param host: elasticsearch host name
        :type host: str
        :param port: elasticsearch port
        :type port: int
        :param api_key: api key to connect to elasticsearch
        :type api_key: str
        :return elasticsearch client object
        """
        try:
            es_client = elasticsearch.Elasticsearch([f"http://{host}:{port}"], api_key=api_key)
        except ElasticsearchError as error:
            logger.error(error)
            raise
        except Exception as error:
            msg = f"Unknown error encountered:  {error.__class__.__name__}: {error}"
            logger.error(msg)
            raise ElasticsearchError(msg)

        return es_client

    @staticmethod
    def es_response_to_json(response: elastic_transport) -> List[dict]:
        """Given an elasticsearch response document, parse the doc for the metrics that are used when
        inserting into elasticsearch.
        """
        jsondata = []
        if "hits" in response:
            for hit in response["hits"]["hits"]:
                jsondata.append(hit["_source"])
        return jsondata

    def query_elasticsearch(self, index: str, query: dict) -> elastic_transport:
        """
        Search Elastic Search (ES) DB

        :param index: elasticsearch index to query
        :type index: str
        :param query: elasticsearch query dictionary
        :type query: dict
        :return: ES search response
        :rtype: elastic_transport object
        """
        try:
            response = self.es_client.search(index=index, **query)
        except elasticsearch.TransportError as error:
            logger.error(error)
            raise ElasticsearchError(error)

        return response

    def insert_into_elasticsearch(self, doc: dict) -> elastic_transport:
        """Given a dict, insert into elasticsearch.

        :param doc: dictionary of the doc to insert (index) into elasticsearch.
        :type doc: dic
        :return elasticsearch response
        """
        index = self.es_index["runs"]
        id = doc.get("run_id", None)

        if not id:
            msg = "ERROR: failed to insert into elasticsearch. JSON doc does not contain run_id."
            logger.error(msg)
            raise ElasticsearchError(msg)

        try:
            response = self.es_client.index(index=index, id=id, document=doc)
        except ElasticsearchError as error:
            logger.error(error)
            raise
        except Exception as error:
            msg = f"ERROR: Unknown error encountered:  {error.__class__.__name__}: {error}"
            logger.error(error)
            raise ElasticsearchError(error)

        return response

    def search_by_run_id(self, run_id: str) -> List[dict]:
        """Given a jaws run id, search elasticsearch for run metadata for that run and return dict of the metadata.

        :param run_id: jaws run id
        :type run_id: int
        :return run metadata dictionary
        """
        index = self.es_index["runs"]
        query = {
            "query": {
                'match': {"run_id": run_id}
            },
            "size": 10000
        }
        response = self.query_elasticsearch(index, query)
        return self.es_response_to_json(response)

    def search_perf_metrics_by_task(self, run_id: str, task_name: str) -> List[dict]:
        """Given a jaws run id and a task name, lookup all performance metric docs in elasticsearch for the run id
        and task name. Return a list of performance metric dictionaries.

        :param run_id: jaws run id
        :type run_id: int
        :param task_name: jaws task name
        :type task_name: str
        :return list of performance metric dicts
        """

        index = self.es_index["performance_metrics"]
        query = {
            'query': {
                'bool': {
                    'must': [
                        {'match': {'jaws_run_id': run_id}},
                        {'match': {'task_name': task_name}},
                    ]
                }
            },
            "sort": [
                {"@timestamp": "asc"}
            ],
            'size': 10000
        }
        response = self.query_elasticsearch(index, query)
        return self.es_response_to_json(response)

    def search_run_by_dates(self, start_date, end_date, size=10000, has_no_perf_metrics_only=False):
        """Given a starting and end date in the format YYYY-MM-DD, query elasticsearch and return all docs
        within that date range. If has_no_perf_metrics_only=True, only return docs that doesn't have performance
        metrics.

        :param start_date: starting date in format YYYY-MM-DD
        :type start_date: str
        :param end_date: ending date in format YYYY-MM-DD
        :type end_date: str
        :param size: max size of docs to return from elasticsearch
        :type size: int
        :param has_no_perf_metrics_only: set to True to return only docs that don't have performance mertrics.
        :type has_no_perf_metrics_only: bool
        :return list of run metadata dicts
        """
        # start_date and end_date must be in format YYYY-MM-DD
        index = self.es_index["runs"]
        date_range = {
            "range": {
                "submitted": {
                    "gte": f"{start_date} 00:00:00",
                    "lte": f"{end_date} 23:59:59"
                }
            }
        }
        query = {
            "query": {},
            "size": size
        }

        if has_no_perf_metrics_only:
            query["query"] = {
                "bool": {
                    "must_not": [
                        {
                            "exists": {
                                "field": "cpu_max"
                            }
                        }
                    ],
                    "filter": [date_range]
                }
            }
        else:
            query["query"] = date_range

        response = self.query_elasticsearch(index, query)
        return self.es_response_to_json(response)
