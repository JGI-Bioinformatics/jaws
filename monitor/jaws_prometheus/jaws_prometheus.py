#!/usr/bin/env python
"""
This program serves as a prometheus custom exporter (otherwise known as a custom collector) for
monitoring jaws components.

The program runs as a daemon and serves as a light-weight http server. The prometheus server
makes an HTTP request to this application to gather the reported metrics and stores the metrics into
the prometheus PromQL database for querying and web display using Grafana.

A config file is used to both define the services to monitor and the mechanism to query the
service for statuses. The supported mechanisms perform either an HTTP or RabbitMQ request to
the service in order to determine it's health or aquire metrics. The statuses are reported
in a specific format designed for the prometheus server to consume.

The config file is either pass in as command line input or set as an environment variable
using JAWS_MONITOR_CONFIG environment variable name.

The structure of the config file is in ini format and looks for the following sections:

[MONITOR]  # required section
port = port for this application

[REST_SERVICES]  # (optional) returns 0 or 1 designating whether the http request failed
                 # or is successful.
name_of_metric = http://somehost:someport
...

[DISK_MONITOR]  # (optional) returns float from http request for returning disk free percent
as float.
name_of_metric = http://somehost:someport
...

[SUPERVISOR_STATUS]  # (optional) returns float from http request for returning supervisord pid
as float.
name_of_metric = http://somehost:someport?config=/path/supervisor.conf&cmd=/path/supervisord
...

[SUPERVISOR_PROCESS]  # (optional) returns float from http request for returning the status of each
process that supervisord is managing.
name_of_metric = http://somehost:someport?config=/path/supervisor.conf&cmd=/path/supervisorctl
...

[RMQ:name_of_metric]   # (optional) performs an rmq/rpc request and returns 0 or 1 designating rmq request
                       # failed or is successful.
                       # The return value must be in JSON_RPC format with the key='result' and
                       # a value of either True or False.
user: rmq_user
password: rmq_pwd
host: rmq_host
port: rmq_port
vhost: rmq_vhost
queue: rmq_queue

...

"""

import os
import requests
import argparse
import time
import json
import configparser
import time
from prometheus_client import start_http_server, Gauge
from pprint import pprint
from jaws_prometheus.config import Configuration
from jaws_rpc import rpc_client

# ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
# import sys
# sys.path.insert(0, os.path.join(ROOT_DIR, '../jaws_prometheus'))
# sys.path.insert(0, os.path.join(ROOT_DIR, '../../rpc'))
# from jaws_rpc import rpc_client
# from config import Configuration



def get_args():
    """Parse command line arguments.

    :param none
    :type none
    :return: argparser object
    :rtype: object
    """

    prog_desc = '''
    '''
    parser = argparse.ArgumentParser(description=prog_desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('config_file', help='jaws monitor config file')
    return parser.parse_args()


def is_http_status_valid(status_code):
    """Check the HTTP status code and return 1 if response is successful or 0 if response
    is bad.

    :param status_code: HTTP status code
    :type status_code: int
    :return 0 or 1
    :rtype: int
    """

    return 1 if status_code >= 200 and status_code < 300 else 0


def set_prometheus_metric(proms, name, value):
    """Given a metric value, report the metric to the prometheus gauge object. If the guage object
    doesn't exist for the given name in the proms dictionary, one is created.

    :param proms: dictionary containing prometheus gauge objects
    :type proms; dictionary
    :param name: name of metric to report
    :type name: string
    :param value: metric value to report
    :type value: integer
    :return: None
    """

    if name not in proms:
        proms[name] = Gauge(name, f"Monitor JAWS {name}")
    print(f"{value}\t{name}")
    proms[name].set(value)


def http_request(url, **rkwargs):
    """Performs an HTTP request and lookup the HTTP return status. If status 200-299, return 1
    else return 0. If kwargs[tofloat] = True, return float form of value from http request.

    :param url: url of HTTP request
    :type url: string
    :return status: http status code
    :rtype status: int
    :return jsondata: json output from http request
    :rtype jsondata: dict
    """

    jsondata = {}
    try:
        if '?' in url:
            url, p = url.split('?')
            rkwargs['params'] = dict(item.split("=") for item in p.split("&"))
        r = requests.get(url=url, **rkwargs)
    except requests.exceptions.RequestException as err:
        return 500, {}
    except Exception as err:
        return 500, {}
    if r.status_code != 200:
        return 500, {}
    try:
        jsondata = r.json()
    except ValueError as err:
        raise ValueError("%s error occurred: %s"%(type(err).__name__, err))

    status = r.status_code

    return status, jsondata


def rpc_request(entries):
    """Create an rpc_client object to connect to the RabbitMQ server, vhost and queue.
    If the connection to the RabbitMQ server fails, return None.

    :param entries: A dictionary containing the rabbitmq connection information. Ex:
    {
        user: RabbitMQ user
        password: RabbitMQ password
        host: RabbitMQ host
        port: RabbitMQ port
        vhost: vhost name
        queue: queue name
    }
    :type entries: dict
    :return: rpc_client object with connection to RabbitMQ
    :rtype rpc_client object
    """

    try:
        rpc = rpc_client.RpcClient(entries)
    except Exception as err:
        rpc = None
    return rpc


def report_rest_services(config, proms):
    """Performs an HTTP request for the given url, checks the http status code for valid ranges and returns 0 or 1,
    with 0 for bad response and 1 for successful response. Also returns the output of the http request.

    :param config: config object
    :type config: jaws_prometheus.config.Configuration object
    :param proms: dictionary containing the prometheus gauge object for reporting metrics
    :type proms; dictionary
    :return: None
    """
    for name, url in config.get_config_section('REST_SERVICES'):
        status_code, _ = http_request(url)
        status = is_http_status_valid(status_code)
        set_prometheus_metric(proms, name, status)


def report_disk_free(config, proms):
    """Performs an HTTP request for the given url. If status is successful, parse json
    output from disk monitor and return disk free percentage as a float. Assumes
    url returns json structure in the format:
    {
        'disk_free_pct': numeric
    }
    If HTTP response is unsuccessful, return -1.

    :param config: config object
    :type config: jaws_prometheus.config.Configuration object
    :param proms: dictionary containing the prometheus gauge object for reporting metrics
    :type proms; dictionary
    :return: None
    """
    for name, url in config.get_config_section('DISK_MONITOR'):
        disk_free_pct = -1
        status_code, jsondata = http_request(url)
        status = is_http_status_valid(status_code)
        if status:
            disk_free_pct = float(jsondata.get('disk_free_pct', 0))
        set_prometheus_metric(proms, name, disk_free_pct)


def report_rmq_services(config, proms):
    """Performs an RPC client request for the give rpc client object from the rpc_client module to
    check for the status of the service or component.
    The input rpc object contains the RMQ connection for a specific channel and queue. The response
    will be reported back by the server process (site, jtm, ...)

    :param config: config object
    :type config: jaws_prometheus.config.Configuration object
    :param proms: dictionary containing the prometheus gauge object for reporting metrics
    :type proms; dictionary
    :return: None
    """

    for name, entries in config.get_config_section('RMQ:', is_partial_name=True):
        status = 0
        rpc = rpc_request(entries)
        if rpc:
            try:
                response = rpc.request("server_status")
            except rpc_client.ConnectionError:
                status = 0
            if response and 'result' in response and response['result'] is True:
                status = 1
            else:
                status = 0
        set_prometheus_metric(proms, name, status)


def report_supervisor_pid(config, proms):
    """Performs a REST call to site-monitor to get the supervisord process id. If alive, site-monitor
    returns a non-zero number in a json structure like:
    {
        'pid': number
    }

    :param config: config object
    :type config: jaws_prometheus.config.Configuration object
    :param proms: dictionary containing the prometheus gauge object for reporting metrics
    :type proms; dictionary
    :return: None
    """

    for name, url in config.get_config_section('SUPERVISOR_STATUS'):
        status = 0
        http_status_code, jsondata = http_request(url)
        if is_http_status_valid(http_status_code):
            status = jsondata.get('pid', 0)
        set_prometheus_metric(proms, name, status)


def report_supervisor_processes(config, proms):
    """Performs a REST call to site-monitor to get the status of all processes managed by the supervisord.
    The site-monitor returns a dictionary where the key is the name of the process, the value is the status
    of the process (0=down, 1=up).

    :param config: config object
    :type config: jaws_prometheus.config.Configuration object
    :param proms: dictionary containing the prometheus gauge object for reporting metrics
    :type proms; dictionary
    :return: None
    """

    for name, url in config.get_config_section('SUPERVISOR_PROCESS'):
        http_status_code, jsondata = http_request(url)
        if not is_http_status_valid(http_status_code):
            continue
        for process in jsondata:
            # Prepend the name of the metric from the config file to the name of the supervisor process.
            metric_name = f"{name}_{process}"

            # prometheus doesn't like dashes in the metric name so we need to convert this
            # to underscores.
            metric_name = metric_name.replace('-', '_')

            status = jsondata[process]
            set_prometheus_metric(proms, metric_name, status)


def main():
    """Main routine to lookup the status of each service defined in the config file and
    report the metrics for the prometheus server to consume.

    :param: none
    :type: none
    :rtype: none
    """

    # Parse args from command-line
    args = get_args()
    config_file = args.config_file
    time_delay = 30 # 30 seconds

    # Parse config file and define services to monitor
    config= Configuration(config_file)

    # Lookup port number for this app from the config file
    port = config.get_port()
    if not port:
        raise SystemExit("Failed to get port number for jaws_prometheus web server")

    # Start http server using the prometheus client module
    start_http_server(port)

    # Here, we define a dictionary to store the prometheus client objects for reporting
    # metrics, one for each metric that we want to capture. Because prometheus only
    # allows the instantiation to be done once, we'll store the client object in this
    # dict as one is created and reuse ones that have already been created. The client
    # object are created when calling the set_prometheus_metric function if not already
    # defined.
    proms = {}

    # Lookup status of each service and report metrics.
    while True:
        report_rest_services(config, proms)
        report_rmq_services(config, proms)
        report_disk_free(config, proms)
        report_supervisor_pid(config, proms)
        report_supervisor_processes(config, proms)

        time.sleep(time_delay)


if __name__ == '__main__':
    main()
