import argparse
import socket
import sys

from jaws_jtm.lib.jtminterface import JtmInterface
from jaws_jtm.lib.config import COMPUTE_RESOURCES


def check_worker():
    desc = u"jtm-check-worker"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-tp", "--task-pool", dest="task_pool")  # if custom pool
    parser.add_argument("-s", "--slurm", dest="slurm_info")
    parser.add_argument(
        "-cl",
        "--cluster",
        dest="cluster_name",
        required=False,
        choices=COMPUTE_RESOURCES,
    )
    parser.add_argument(
        "-l",
        "--loglevel",
        help="Set loglevel (default=info).",
        dest="log_level",
        required=False,
        default=None,
    )
    args = parser.parse_args()

    add_info = socket.gethostname().replace(".", "_")
    if args.cluster_name:
        add_info = args.cluster_name

    if args.slurm_info:  # todo: return slurm info
        ret = JtmInterface("check_worker", info_tag=add_info).call(
            task_pool=args.task_pool,
            slurm_info=args.task_pool,  # todo, not used.
            jtm_host_name=args.cluster_name,
            log_level=args.log_level,
        )

    else:
        ret = int(
            JtmInterface("check_worker", info_tag=add_info).call(
                task_pool=args.task_pool,
                jtm_host_name=args.cluster_name,
                log_level=args.log_level,
            )
        )
        print(ret)
    sys.exit(0) if ret > 0 else sys.exit(1)
