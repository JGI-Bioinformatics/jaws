#! /usr/bin/env python
#
# Check the number of workers
# exit >0 if found
#       0 if no worker found
#
import argparse
import sys

from jaws_jtm.lib.jtminterface import JtmInterface
from jaws_jtm.config import COMPUTE_RESOURCES


def remove_pool():
    desc = u"jtm-remove-pool"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-p", "--pool_name",
                        dest="task_pool",
                        required=True)  # if custom pool
    parser.add_argument("-cl", "--cluster",
                        dest="cluster_name",
                        required=False,
                        choices=COMPUTE_RESOURCES)
    parser.add_argument("-l", "--loglevel",
                        help="Set loglevel (default=info).",
                        dest="log_level", required=False,
                        default=None)
    args = parser.parse_args()
    ret = int(JtmInterface('remove_pool', info_tag=args.task_pool).call(task_pool=args.task_pool,
                                                                        jtm_host_name=args.cluster_name,
                                                                        log_level=args.log_level))
    print("removed" if ret == 1 else "failed")
    sys.exit(0) if ret > 0 else sys.exit(1)
