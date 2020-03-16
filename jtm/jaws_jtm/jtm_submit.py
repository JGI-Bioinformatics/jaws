#! /usr/bin/env python
#
# JtmInterface retunrs code
#   'task_id' if successfully queued
#
# jtm-submit exits with code 0 if successfully submitted
#                            1 if not
#
# USAGE
# jtm-submit -f <json_file>
# jtm-submit -j '<json_str>'
#
# to send task to other cluster
# jtm-submit -f <json_file> -cl <cluster_name>
#
import argparse
import sys

from jaws_jtm.lib.jtminterface import JtmInterface
from jaws_jtm.config import COMPUTE_RESOURCES, \
    QOS_LIST, \
    CORI_QOS, \
    NWORKERS, \
    NNODES, \
    MEMPERNODE, \
    CORI_CONSTRAINT, \
    CORI_CHARGE_ACCNT
from jaws_jtm.lib.run import eprint


def submit():
    desc = u"jtm-submit"
    parser = argparse.ArgumentParser(description=desc)
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("-cl", "--cluster",
                        dest="cluster_name",
                        choices=COMPUTE_RESOURCES)
    group.add_argument("-cr", "--cromwell",
                       dest="cromwell")  # generic for command submission ex) jtm-submit -cr "ls"
    group.add_argument("-f", "--task-json-file",
                       dest="taskFile")
    group.add_argument("-j", "--task-json-str",
                       dest="task_in_json_str")
    parser.add_argument("-l", "--loglevel",
                        help="Set loglevel (default=info).",
                        dest="log_level",
                        required=False,
                        default=None)

    # resource requirement
    # This will be used for create a separate pool of jtm-worker(s)
    parser.add_argument("-A", "--account",
                        dest="account",
                        default=CORI_CHARGE_ACCNT)
    parser.add_argument("-c", "--cpu",
                        help="Number of cores",
                        type=int,
                        dest="num_cpus")
    parser.add_argument("-C", "--constraint",  # only for nersc
                        help="Set the architecture to Haswell or KNL on Cori",
                        dest="constraint",
                        choices=["haswell", "knl", "skylake"],
                        default=CORI_CONSTRAINT)
    parser.add_argument("-jid", "--job-id",
                        help="Unique Cromwell job id with step name",
                        dest="cromwell_job_name",
                        required=False)
    parser.add_argument("-m", "--memory",
                        help="Memory request per node",
                        dest="memory",
                        default=MEMPERNODE)
    parser.add_argument("-N", "--num-nodes",
                        help="Number of nodes for the pool. Default=1.",
                        dest="num_nodes",
                        type=int,
                        default=NNODES)
    parser.add_argument("-nwpn", "--num-workers-per-node",
                        help="Number of worker per node. Default=1",
                        dest="num_worker_per_node",
                        required=False,
                        type=int,
                        default=NWORKERS)
    parser.add_argument("-p", "--pool-name",
                        dest="pool_name",
                        required=True)
    parser.add_argument("-q", "--qos",
                        dest="qos",
                        choices=QOS_LIST,
                        default=CORI_QOS)
    parser.add_argument("-s", "--shared",
                        help="Shared workers.",
                        dest="shared_worker",
                        type=int,
                        default=1)
    parser.add_argument("-t", "--time",
                        dest="job_time")

    args = parser.parse_args()

    # if any(vars(args).values()):
    #     args = parser.parse_args(['-h'])

    if args.job_time:
        assert args.num_cpus and args.memory and \
               args.num_nodes and args.pool_name, \
               "USAGE: runtime (-t) should be set with 'num_cpus' (-c), memory (-m), node (-nn), and pool name (-p)."

    if args.pool_name == "default":
        args.pool_name = "small"

    json_task_str = ""
    if args.cromwell:
        json_task_str = args.cromwell
        # assert args.job_time and args.cpus and args.memory and args.pool_name
    elif args.task_in_json_str:
        json_task_str = args.task_in_json_str

    add_info = None
    if args.pool_name:
        add_info = args.pool_name

    ret = int(JtmInterface('task', info_tag=add_info).call(task_file=args.taskFile,
                                                           task_json=json_task_str,
                                                           task_id=0,
                                                           jtm_host_name=args.cluster_name,
                                                           job_time=args.job_time,
                                                           node_mem=args.memory,
                                                           num_core=args.num_cpus,
                                                           pool_name=args.pool_name,
                                                           log_level=args.log_level,
                                                           shared=args.shared_worker,
                                                           nwpn=args.num_worker_per_node,
                                                           node=args.num_nodes,
                                                           job_id=args.cromwell_job_name,
                                                           constraint=args.constraint,
                                                           qos=args.qos,
                                                           account=args.account
                                                           ))

    if ret == -5:
        eprint("jtm-submit: invalid task or runtime definition.")
        return 1
    elif ret == -88:
        eprint("jtm-submit: command timeout.")
        return -1

    print("JTM task ID %d" % ret)
    return 0 if ret != 0 else 1


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    sys.exit(submit())
