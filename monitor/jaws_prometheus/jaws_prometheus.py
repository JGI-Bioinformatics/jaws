#!/usr/bin/env python
"""
This program serves as a prometheus custom exporter (otherwise known as a custom collector) for
monitoring jaws components.

The program runs as a daemon and serves as a light-weight web server. The prometheus server
makes a HTTP request to this application to gather the reported metrics and stores the metrics into
the prometheus PromQL database for querying and web display using Grafana.

A config file is used to both define the services to monitor and the mechanism to query the
service for statuses. The supported mechanisms perform either an HTTP or RabbitMQ request to
the service in order to determine it's health. The statuses are reported as metrics in
a format for the prometheus server to consume.

The config file is either pass in as command line input or set as an environment variable
using JAWS_MONITOR_CONFIG environment variable name.

The structure of the config file is in ini format and looks for the following sections:

[MONITOR]    # required section
port = port for this application

[REST:BOOL]  # (optional) returns 0 or 1 designating whether the http request failed
             # or is successful.
name_of_service_1 = http://somehost:someport
name_of_service_2 = http://somehost:someport
...

[REST:FLOAT] # (optional) returns float from http request for returning disk free percent.
name_of_service_3 = http://somehost:someport
name_of_service_4 = http://somehost:someport
...

[RMQ:name_of_service_5]  # (optional) performs an rmq/rpc request and returns 0 or 1 designating rmq request
                         # failed or is successful.
                         # The return value must be in JSON_RPC format with the key='result' and
                         # a value of either True or False.
user: rmq_user
password: rmq_pwd
host: rmq_host
port: rmq_port
vhost: rmq_vhost
queue: rmq_queue

[RMQ:name_of_service_6] # (optional)
...

"""

import os
import requests
import argparse
import time
import configparser
from prometheus_client import start_http_server, Gauge

# ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
# sys.path.insert(0, os.path.join(ROOT_DIR, '../../rpc'))
from jaws_rpc import rpc_client


class Config():
    """Class to read and get config file entries.
    """
    def __init__(self, config_file):
        """Read the given config file and create a config parser object containing the
        contents of the config file.

        :param config_file: file path and name of config file
        :type config_file: str
        :return: none
        """
        if not config_file:
            config_file = os.environ.get('JAWS_MONITOR_CONFIG')
        config = configparser.ConfigParser()
        try:
            config.read(config_file)
        except Exception:
            raise
        self.config = config

    def get_services(self):
        """Parse the config entries and create a dictionary containing the following keys
        for each component to be monitored:

        function=name of function to call to get status
        params=parameter(s) used to pass into function for getting status

        The functions here are either for http requests or for rabbitmq requests.

        Example of returned dict:
        services = {
            "central_prod": {
                "function": http_respones
                "params": "http://jaws.lbl.gov:5003/api/v2/status"
            }
            "site_cori_prod": {
                "function":rmq_request
                "params": rpc_client.RpcClient()
            }
        }

        :param none
        :type none
        :return: none
        """

        services = {}

        # Add monitoring for services that uses RMQ/RPC
        for section in self.config.sections():
            # Add monitoring for services that uses REST
            if section == 'REST:BOOL':
                for name in self.config[section]:
                    services[name] = {
                        'function': http_response,
                        'params': self.config['REST:BOOL'].get(name)
                    }

            # Add monitoring for disk quota using REST
            elif section == 'REST:FLOAT':
                for name in self.config[section]:
                    services[name] = {
                        'function': disk_free,
                        'params': self.config['REST:FLOAT'].get(name)
                    }

            # Add monitoring for services that uses RMQ/RPC
            elif section.startswith('RMQ:'):
                services[section.lower()] = {
                    'function': rmq_request,
                    'params': rpc_client.RpcClient(self.config[section]),
                }

        return services

    def get_port(self):
        """Parse the config entries to get the port number for serving up this web server using
        the prometheus client module.

        :param none
        :type none
        :return: none
        """
        port = self.config['MONITOR'].get('port')
        if port:
            return int(port)
        else:
            return None


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


def rmq_request(rpc):
    """Perform an RPC client request for the give rpc client object from the rpc_client module to
    check for the status of the service or component.
    The input rpc object contains the RMQ connection for a specific channel and queue. The response
    will be reported back by the server process (site, jtm, ...)

    :param rpc: rpc_client object for rabbitMQ requests
    :type rpc: object
    :return: 0 or 1
    :rtype: int
    """
    try:
        response = rpc.request("server_status")
    except rpc_client.ConnectionError:
        return 0
    if response and 'result' in response and response['result'] is True:
        return 1
    else:
        return 0


def create_proms(services):
    """For each service or component to be monitored (defined in the services dict config entries),
    create a prometheus gauge object to report metrics to the prometheus server.

    :param services: dictionary defining the services to monitor.
    :type services: dict
    :return: dictionary where key=name of service, value=prometheus gauge object.
    :rtype: object
    """
    # Register and create prometheus objects for each service to store metrics.
    proms = {}
    for name in services:
        proms[name] = Gauge(name, f"Monitor JAWS {name}")
    return proms


def disk_free(*args, **rkwargs):
    """Perform an HTTP request looking for a float like string as a return value. If the return
    is successful, convert the string format to a float and return the float value.

    :param *args: positional arguments; 1st element must be a url to the server program that
        reports the disk free percentage
    :type *args: list
    :return retval: disk free percentage or -1 if http request fails.
    :rtype retval: float | int
    """
    if not len(args) == 1:
        raise SystemExit("Missing url parameter when calling rest_request()")
    url = args[0]

    retval = 0

    try:
        r = requests.get(url=url, **rkwargs)
    except Exception:
        return -1

    if r.status_code >= 200 and r.status_code < 300:
        try:
            retval = float(r.text)
        except ValueError:
            retval = -1
    else:
        retval = -1

    return retval


def http_response(*args, **rkwargs):
    """Performs an HTTP request and lookup the HTTP return status. If status 200-299, return 1
    else return 0

    :param *args: positional arguments; 1st element must be a url to the server program that
        reports the disk free percentage
    :type *args: list
    :return retval: 0 or 1
    :rtype retval: int
    """
    if not len(args) == 1:
        raise SystemExit("Missing url parameter when calling http_response()")
    url = args[0]

    try:
        r = requests.get(url=url, **rkwargs)
    except Exception:
        return 0
    return 1 if r.status_code >= 200 and r.status_code < 300 else 0


def process_request(services, proms, t=60):
    """Lookup the status for each service defined in the services dict (derived from the config
    file entries) and report the metrics using the prometheus client module.

    :param services: dictionary defining the services to monitor.
    :type services: dict
    :param proms: dictionary where key=name of service, value=prometheus gauge object
    :type proms: dict
    :param t: sleep time in seconds
    :type t: int
    :rtype: none
    """

    for name in services:
        func = services[name].get('function')
        if not func:
            raise SystemExit("Missing func key in services dict")

        params = services[name].get('params')
        status = func(params)
        proms[name].set(status)
        print(name, params, status)

    time.sleep(t)


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

    # Parse config file and define services to monitor
    config = Config(config_file)
    services = config.get_services()

    port = config.get_port()
    if not port:
        raise SystemExit("Failed to get port number for jaws_prometheus web server")

    # Start http server using the prometheus client module
    start_http_server(port)

    # Create prometheus gauge objects for each service to monitor.
    proms = create_proms(services)

    # Lookup status of each service and report status.
    while True:
        process_request(services, proms)


if __name__ == '__main__':
    main()
