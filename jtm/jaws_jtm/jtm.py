#!/usr/bin/env python
import sys
import os
import click
import socket
import json
import csv
import time
import re

from jaws_jtm.config import JtmConfig
from jaws_jtm.jtm_manager import manager as jtmmanager
from jaws_jtm.jtm_worker import worker as jtmworker
from jaws_jtm.lib.jtminterface import JtmInterface
from jaws_jtm.common import logger
from jaws_jtm.lib.run import run_sh_command


class Mutex(click.Option):
    def __init__(self, *args, **kwargs):
        self.not_required_if: list = kwargs.pop("not_required_if")
        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs["help"] = (
            kwargs.get("help", "")
            + "Option is mutually exclusive with "
            + ", ".join(self.not_required_if)
            + "."
        ).strip()
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt: bool = self.name in opts
        for mutex_opt in self.not_required_if:
            if mutex_opt in opts:
                if current_opt:
                    raise click.UsageError(
                        "Illegal usage: '"
                        + str(self.name)
                        + "' is mutually exclusive with "
                        + str(mutex_opt)
                        + "."
                    )
                else:
                    self.prompt = None
        return super(Mutex, self).handle_parse_result(ctx, opts, args)


@click.group()
@click.option("--debug", is_flag=True, default=False)
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.pass_context
def cli(ctx: object, debug: bool, config_file: str):
    # click.echo("Debug mode is %s" % ("on" if debug else "off"))
    if config_file:
        config = JtmConfig(config_file=config_file)
    else:
        config = JtmConfig()
    ctx.obj = {"config_file": config_file, "config": config, "debug": debug}


@cli.command()
@click.option("-ld", "--log_dir", help="Custom log directory", required=False)
@click.option(
    "-r",
    "--show_resource_log",
    help="Show resource usage log from worker(s)",
    default=False,
    is_flag=True,
)
@click.pass_context
def manager(ctx: object, log_dir: str, show_resource_log: bool) -> int:
    """
    JTM Manager Click wrapper

    :param ctx:
    :param log_dir: custom log dir
    :param show_resource_log: show/not show resource usage log
    :return:
    """
    sys.exit(jtmmanager(ctx, log_dir, show_resource_log))


@cli.command()
@click.option(
    "-hb",
    "--heartbeat_interval",
    help="Heartbeat interval in second (default=10).",
    type=int,
)
@click.option("-jd", "--job_script_dir_name", help="Slurm batch job file path.")
@click.option("-ld", "--log_dir", help="Custom log directory", required=False)
@click.option(
    "-p",
    "--pool_name",
    help="Set user-defined pool name. This should be same with task's 'pool' name.",
    required=True,
)
@click.option(
    "-dr",
    "--dry_run",
    help="Dry run. To print batch job script on screen",
    default=False,
    is_flag=True,
)
@click.option("-j", "--slurm_job_id", help="Slurm job ID", type=int, default=0)
@click.option(
    "-wt",
    "--worker_type",
    help="Worker type = [manual | dynamic]",
    default="manual",
    type=click.Choice(["manual", "dynamic"], case_sensitive=False),
)
@click.option("-cl", "--cluster", help="Site name")
@click.option(
    "-ctr", "--clone_time_rate", help="Cloning time rate (OBSOLETE)", type=float
)
@click.option(
    "-nwpn",
    "--num_worker_per_node",
    help="Set number of workers per a node",
    required=False,
    type=int,
    default=1,
)
@click.option("-wi", "--worker_id", help="Unique worker ID")
@click.option("-A", "--charging_account", help="Slurm charging account.")
@click.option("-N", "--nnodes", help="Slurm number of nodes", type=int)
@click.option("-c", "--cpus_per_task", help="Slurm number of cpus", type=int)
@click.option(
    "-C", "--constraint", help="Slurm constraint (default: haswell)", required=False
)
@click.option("-m", "--mem", help="Slurm real memory required per node")
@click.option(
    "-mc", "--mem_per_cpu", help="Slurm minimum memory required per allocated CPU"
)
@click.option("-P", "--partition", help="Slurm partition name")
@click.option("-q", "--qos", help="Slurm quality of service")
@click.option(
    "-r",
    "--show_resource_log",
    help="Show resource usage log from worker(s)",
    default=False,
    is_flag=True,
)
@click.option("-t", "--job_time", help="Slurm Job time (hh:mm:ss)")
@click.pass_context
def worker(
    ctx: object,
    heartbeat_interval: int,
    log_dir: str,
    job_script_dir_name: str,
    pool_name: str,
    dry_run: bool,
    slurm_job_id: int,
    worker_type: str,
    cluster: str,
    clone_time_rate: float,
    num_worker_per_node: int,
    worker_id: str,
    charging_account: str,
    nnodes: int,
    cpus_per_task: int,
    constraint: str,
    mem: str,
    mem_per_cpu: str,
    qos: str,
    partition: str,
    job_time: str,
    show_resource_log: bool,
) -> int:
    """
    JTM Worker Click wrapper

    :param ctx:
    :param heartbeat_interval: worker HB interval to the manager
    :param log_dir: custom log directory
    :param job_script_dir_name: custom sbatch script saving directory
    :param pool_name: worker pool name
    :param dry_run:
    :param slurm_job_id: SLURM job id
    :param worker_type: manual or dynamic
    :param cluster: destination cluster name
    :param clone_time_rate: OBSOLETE
    :param num_worker_per_node:
    :param worker_id:
    :param charging_account: SLURM charging account
    :param nnodes: number of nodes
    :param cpus_per_task: number of cpus
    :param constraint: SLURM constraint
    :param mem: memory request
    :param mem_per_cpu: memory request per cpu
    :param qos:
    :param partition:
    :param job_time: wallclocktime
    :return:
    """
    config = ctx.obj["config"]
    if cluster and cluster.upper() not in config.constants.COMPUTE_RESOURCES:
        logger.error(f"Unsupported site name specified: {cluster}")
        sys.exit(1)

    sys.exit(
        jtmworker(
            ctx,
            heartbeat_interval,
            log_dir,
            job_script_dir_name,
            pool_name,
            dry_run,
            slurm_job_id,
            worker_type,
            cluster,
            clone_time_rate,
            num_worker_per_node,
            worker_id,
            charging_account,
            nnodes,
            cpus_per_task,
            constraint,
            mem,
            mem_per_cpu,
            qos,
            partition,
            job_time,
            show_resource_log,
        )
    )


@cli.command()
@click.option("-cl", "--cluster", help="Site name")
@click.option(
    "-cmd",
    "--command",
    help="User task command",
    cls=Mutex,
    not_required_if=["task_file"],
)
@click.option(
    "-f",
    "--task_file",
    help="User task in a json file",
    cls=Mutex,
    not_required_if=["command"],
)
@click.option(
    "-A",
    "--account",
    help="Slurm charging account name",
)
@click.option("-c", "--ncpu", help="Slurm number of cores", type=int)
@click.option("-C", "--constraint", help="Set the architecture for SLURM request")
@click.option(
    "-jid",
    "--cromwell_job_id",
    help="Unique Cromwell job id with step name",
    required=False,
)
@click.option("-m", "--memory", help="Slurm memory request per node", type=float)
@click.option(
    "-N",
    "--nnodes",
    help="Slurm number of nodes for the pool. Default=1.",
    type=int,
    default=1,
)
@click.option(
    "-nwpn",
    "--num_worker_per_node",
    help="Number of worker per node. Default=1",
    required=False,
    type=int,
    default=1,
)
@click.option("-p", "--pool_name", help="User pool name", required=True)
@click.option("-P", "--partition", help="Slurm partition name")
@click.option("-q", "--qos", help="Set the QOS for SLURM request")
@click.option("-s", "--shared", help="Shared/non-shared worker", type=int, default=1)
@click.option(
    "-t",
    "--job_time",
    help="Job time (hh:mm:ss)",
)
@click.pass_context
def submit(
    ctx: object,
    task_file: str,
    cluster: str,
    command: str,
    pool_name: str,
    account: str,
    ncpu: int,
    constraint: str,
    num_worker_per_node: int,
    cromwell_job_id: str,
    memory: float,
    nnodes: int,
    qos: str,
    partition: str,
    shared: int,
    job_time: str,
) -> int:
    """
    JtmInterface returns 'task_id' if successfully queued
    jtm submit exits with code 0 if successfully submitted
                               1 if not

    USAGE
    $ jtm submit -cmd "ls" -p test_pool
    $ jtm submit -f <json_file>

    :param ctx:
    :param task_file:
    :param cluster:
    :param command:
    :param pool_name:
    :param account:
    :param ncpu:
    :param constraint:
    :param num_worker_per_node:
    :param cromwell_job_id:
    :param memory:
    :param nnodes:
    :param qos:
    :param partition:
    :param shared:
    :param job_time:
    :return: exit with 0 if succeeded, 1 if failed
    """
    config = ctx.obj["config"]
    ret_code = 0
    if ctx.obj["debug"]:
        click.echo("Debug mode")

    if shared is None:
        logger.error("jtm submit: invalid task or runtime definition. ")
        logger.error("'shared' option cannot be None")
        sys.exit(1)

    cluster = cluster.lower()
    if cluster is None or cluster.upper() not in config.constants.COMPUTE_RESOURCES:
        logger.error(f"Unsupported site name specified: {cluster}")
        logger.error(
            f"Only {config.constants.COMPUTE_RESOURCES} are supported as site name"
        )
        sys.exit(1)

    def is_time_format(input: str):
        # Allowed format: 12:12:12 02:59:59 77:00:00 12312312:00:00
        patt = re.compile(r"\d+:[0-5][0-9]:[0-5][0-9]")
        if patt.match(input):
            return True
        else:
            return False

    if time is None or not is_time_format(job_time):
        logger.error("jtm submit: invalid task or runtime definition: ")
        logger.error(
            f"The job time is not in an expected format (%H:%M:%S): {job_time}"
        )
        sys.exit(1)

    node_mem_val = int(memory+0.5)
    node_mem = f"{node_mem_val}G"

    add_info = pool_name
    ret = JtmInterface("task", ctx, info_tag=add_info).call(
        task_file=task_file,
        task_json=command,
        task_id=0,
        jtm_host_name=cluster,
        job_time=job_time,
        node_mem=node_mem,
        num_core=ncpu,
        pool_name=pool_name,
        shared=shared,
        nwpn=num_worker_per_node,
        node=nnodes,
        job_id=cromwell_job_id,
        constraint=constraint,
        qos=qos,
        partition=partition,
        account=account,
    )

    if ret is None or ret == -1:
        logger.error("jtm submit: invalid task or runtime definition.")
        ret_code = 1
    else:
        try:
            # ret is str type so convert it
            ret = int(ret)
        except TypeError as e:
            logger.exception(f"Unexpected return value from the manager: {ret} {e}")
            sys.exit(-1)

    if ret == -2:
        logger.error("jtm submit: invalid task.")
        ret_code = 1
    elif ret == -3:
        logger.error(
            "jtm submit: failed to submit a task. No task ID returned from the manager."
        )
        ret_code = 1
    elif ret == -4:
        logger.error(
            "jtm non-task command: failed to get a response from the manager for non-task commands"
        )
        ret_code = 1
    elif ret == -5:
        logger.error("jtm submit: invalid task or runtime definition.")
        ret_code = 1
    elif ret == -88:
        logger.error("jtm submit: jtm command timeout.")
        ret_code = -1

    # IMPORTANT: this line will regexed by Cromwell to extract a valid task ID
    click.echo("JTM task ID %s" % ret)
    sys.exit(ret_code)


@cli.command()
@click.option("-tid", "--task_id", help="JTM task ID", type=int, required=True)
@click.pass_context
def kill(ctx: object, task_id: int) -> int:
    """
    # JtmInterface returns code
    #   0: terminated successfully
    #  -1: no task found
    #   4: the task has already been completed
    #
    # jtm kill exits with 0 if it's successfully terminated.#
    #                     1 else
    #

    :param ctx:
    :param task_id:
    :return: exit with 0 if succeeded, 1 if failed
    """
    try:
        ret = int(JtmInterface("kill", ctx, info_tag=task_id).call_aux(task_id=task_id))
    except TypeError as e:
        logger.exception(
            f"jtm kill: Unexpected return value from the manager: {ret} {e}"
        )
        sys.exit(1)

    if ret == -88:
        logger.error("jtm kill: command timeout.")
        sys.exit(-1)
    elif ret == -5:
        logger.error("jtm kill: task id not found.")
        sys.exit(-1)
    sys.exit(0) if ret == 0 else sys.exit(1)


@cli.command()
@click.option("-tid", "--task_id", help="JTM task ID", type=int, required=True)
@click.pass_context
def isalive(ctx: object, task_id: int) -> int:
    """
    JtmInterface returns
      0 if ready
      1 if queued
      2 if pending
      3 if running
      4 if successfully done
      -1, -2, -3, -4: failed

    jtm isalive exits with 0 if it's in ['ready', 'queued', 'pending', 'running'] status
                          1 if it's done successfully or failed

    click.echo("isalive")

    :param ctx:
    :param task_id:
    :return: exit with 0 if alive, 1 if not
    """
    try:
        ret = int(
            JtmInterface("status", ctx, info_tag=task_id).call_aux(task_id=task_id)
        )
    except TypeError as e:
        logger.exception(
            f"jtm isalive: Unexpected return value from the manager: {ret} {e}"
        )
        sys.exit(1)

    if ret == -88:
        logger.error("jtm isalive: command timeout.")
        sys.exit(-1)
    elif ret in [0, 1, 2, 3]:
        click.echo("yes")
        sys.exit(0)
    else:
        click.echo("no")
        sys.exit(1)


@cli.command()
@click.option("-tid", "--task_id", help="JTM task ID", type=int, required=True)
@click.pass_context
def status(ctx: object, task_id: int) -> int:
    """
    JtmInterface returns
      0 if ready
      1 if queued
      2 if pending
      3 if running
      4 if successfully done
      -1, -2, -3, -4: failed

    jtm status exits with 0 if it's in ['ready', 'queued', 'pending', 'running'] status
                          1 if it's done successfully or failed

    :param ctx:
    :param task_id:
    :return:
    """
    try:
        ret = int(
            JtmInterface("status", ctx, info_tag=task_id).call_aux(task_id=task_id)
        )
    except TypeError as e:
        logger.exception(
            f"jtm status: Unexpected return value from the manager: {ret} {e}"
        )
        sys.exit(1)

    if ret == -88:
        logger.error("jtm status: command timeout.")
        sys.exit(-1)
    reversed_task_status = dict(
        map(reversed, ctx.obj["config"].constants.TASK_STATUS.items())
    )
    click.echo(reversed_task_status[ret])
    sys.exit(0) if ret in [0, 1, 2, 3] else sys.exit(1)


@cli.command()
@click.option("-p", "--pool_name", help="User worker pool name", required=True)
@click.option("-cl", "--cluster", help="Site name")
@click.pass_context
def remove_pool(ctx: object, pool_name: str, cluster: str) -> int:
    """
    Remove pool of workers from HPC
    This actually looks up SLURM IDs for the pool
    and executes scancel

    # Check the number of workers
    # exit >0 if found
    #       0 if no worker found
    #

    :param ctx:
    :param pool_name:
    :param cluster:
    :return: exit with 0 if succeeded, 1 if failed
    """
    config = ctx.obj["config"]
    if cluster and cluster.upper() not in config.constants.COMPUTE_RESOURCES:
        logger.error(f"Unsupported site name specified: {cluster}")
        sys.exit(1)

    try:
        ret = int(
            JtmInterface("remove_pool", ctx, info_tag=pool_name).call_aux(
                task_pool=pool_name, jtm_host_name=cluster
            )
        )
    except TypeError as e:
        logger.exception(f"Unexpected return value from the manager: {ret} {e}")
        sys.exit(1)

    click.echo("removed" if ret == 1 else "failed")
    sys.exit(0) if ret > 0 else sys.exit(1)


@cli.command()
@click.option("-p", "--pool_name", help="User worker pool name", required=True)
@click.option(
    "-s",
    "--slurm_info",
    help="Show SLURM Job Id information",
    default=False,
    show_default=True,
)
@click.option("-cl", "--cluster", help="Site name")
@click.pass_context
def check_worker(ctx: object, pool_name: str, slurm_info: bool, cluster) -> int:
    """
    total number of workers of the site
    if pool_name is specified,
    returns the number of workers for the pool

    :param ctx:
    :param pool_name:
    :param slurm_info:
    :param cluster:
    :return: number of workers in the pool specified
    """
    config = ctx.obj["config"]
    if cluster and cluster.upper() not in config.constants.COMPUTE_RESOURCES:
        logger.error(f"Unsupported site name specified: {cluster}")
        sys.exit(1)

    add_info = socket.gethostname().replace(".", "_")
    if cluster:
        add_info = cluster

    if slurm_info:  # todo: return slurm job id info
        ret = JtmInterface("check_worker", ctx, info_tag=add_info).call_aux(
            task_pool=pool_name,
            slurm_info=slurm_info,  # todo, not used.
            jtm_host_name=cluster,
        )

    else:
        try:
            ret = int(
                JtmInterface("check_worker", ctx, info_tag=add_info).call_aux(
                    task_pool=pool_name, jtm_host_name=cluster
                )
            )
        except TypeError as e:
            logger.exception(
                f"jtm check_worker: Unexpected return value from the manager: {ret} {e}"
            )
            sys.exit(1)

        click.echo(ret)
    sys.exit(0) if ret > 0 else sys.exit(1)


@cli.command()
@click.option("-cl", "--cluster", help="Site name")
@click.pass_context
def check_manager(ctx: object, cluster: str) -> int:
    """
    check if a JTM manager is running on a site

    :param ctx:
    :param cluster:
    :return: exit with 0 if a manager found, 1 if not
    """
    config = ctx.obj["config"]
    mode = config.configparser.get("JTM", "run_mode")

    if cluster and cluster.upper() not in config.constants.COMPUTE_RESOURCES:
        logger.error(f"Unsupported site name specified: {cluster}")
        sys.exit(1)

    ps_cmd = "ps -aef | grep -v grep | grep -v check-manager | grep jtm | grep manager "
    if mode != "test":
        ps_cmd += f"| grep jaws-{mode} "
    else:
        ps_cmd += "| grep test "
    ps_cmd += "| wc -l"
    num_total_procs = 0
    try:
        so, _, ec = run_sh_command(ps_cmd, log=logger, show_stdout=False)
        num_total_procs = int(so.rstrip())
    except TypeError as te:
        logger.exception(te)
        logger.error(so)
        sys.exit(-2)
    except Exception as e:
        logger.exception(e)
        logger.error(ps_cmd)
        sys.exit(-1)
    else:
        # THIS IS THE NUM OF JTM MANAGER PROCS = 7
        if num_total_procs != config.constants.JTM_NUM_PROCS:
            sys.exit(1)

    add_info = socket.gethostname().replace(".", "_")
    if cluster:
        add_info = cluster
    try:
        ret = int(
            JtmInterface("check_manager", ctx, info_tag=add_info).call_aux(
                jtm_host_name=cluster
            )
        )
    except TypeError as e:
        logger.exception(
            f"jtm check_manager: Unexpected return value from the manager: {ret} {e}"
        )
        sys.exit(1)

    if ret is None or ret == -88:
        sys.exit(-1)

    sys.exit(0) if ret != 0 else sys.exit(1)


@cli.command()
@click.option("-tid", "--task_id", help="JTM task ID", type=int, required=True)
@click.pass_context
def resource_log(ctx: object, task_id: int) -> int:
    """
    Find resource log file from LOG_DIR for a task id
    and return a json format if exist

    :param ctx:
    :param task_id:
    :return:  JSON string
    """
    resource_log_file = JtmInterface("resource", ctx, info_tag=task_id).call_aux(
        task_id=task_id
    )
    if isinstance(resource_log_file, int) and resource_log_file == -88:
        logger.error("jtm resource-log: command timeout.")
        sys.exit(-1)
    elif resource_log_file is None or resource_log_file == "None":
        logger.error("jtm resource-log: invalid resource log file name")
        sys.exit(-1)

    # http://www.andymboyle.com/2011/11/02/quick-csv-to-json-parser-in-python/
    if os.path.isfile(resource_log_file):
        f = open(resource_log_file, "rU")
        reader = csv.DictReader(
            f,
            fieldnames=(
                "child_pid",  # 1
                "clone_time_rate",  # 2
                "cpu_load",  # 3
                "end_date",  # 4
                "host_name",  # 5
                "ip_address",  # 6
                "job_time",  # 7
                "life_left",  # 8
                "mem_per_core",  # 9
                "mem_per_node",  # 10
                "num_cores",  # 11
                "num_tasks",  # 12
                "num_workers_on_node",  # 13
                "perc_mem_used",  # 14
                "pool_name",  # 15
                "ret_msg",  # 16
                "rmem_usage",  # 17
                "root_pid",  # 18
                "run_time",  # 19
                "slurm_jobid",  # 20
                "task_id",  # 21
                "vmem_usage",  # 22
                "worker_id",  # 23
                "worker_type",  # 24
                "jtm_host_name",  # 25
                "nwpn",  # 26
            ),
        )
        click.echo(
            """{ "task_id": %d, "resource_log": %s }"""
            % (task_id, json.dumps([row for row in reader]))
        )
    else:
        logger.error("Resource file, %s, not found." % (resource_log_file))
        resource_log_file = None

    sys.exit(0) if resource_log_file is not None else sys.exit(1)


def jtm():
    cli()
