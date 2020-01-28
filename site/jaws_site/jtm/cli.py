import click
from .Config import *
from .JtmInterface import *

@click.group()
def jtm():
    """
    JTM CLI
    """
    pass

@jtm.command()
@click.argument('task_id')
def status(task_id):
    """
    jtm-status exits with 0 if it's in ['ready', 'queued', 'running'] status
        1 if it's done successfully or failed
    """

    assert len(sys.argv) == 2, "USAGE: jtm-status <task_id>"
    ret = int(JtmInterface('status').call(task_id=int(sys.argv[1])))
    if ret == -88:
        eprint("jtm-status: command timeout.")
        sys.exit(-1)
    reversed_TASK_STATUS = dict(list(map(reversed, list(TASK_STATUS.items()))))
    print(reversed_TASK_STATUS[ret])
    sys.exit(0) if ret in [0, 1, 2] else sys.exit(1)


@jtm.command()
@click.option('--task-pool', default=None, help="task pool")
@click.option('--cluster', type=click.Choice(COMPUTE_RESOURCES), help="cluster")
@click.option('--loglevel', default=None, help="Set loglevel (default=info)")
def check_worker(task_pool, cluster, loglevel):
    ret = int(JtmInterface('check_worker').call(task_pool=task_pool,
                                                jtm_host_name=cluster,
                                                log_level=loglevel))
    print(ret)
    sys.exit(0) if ret > 0 else sys.exit(1)


@jtm.command()
@click.option('--cluster', type=click.Choice(COMPUTE_RESOURCES), help="cluster")
def check_manager(cluster):
    ret = int(JtmInterface('check_manager').call(jtm_host_name=cluster))
    if ret == -88:
        sys.exit(-1)
    sys.exit(0) if ret != 0 else sys.exit(1)

@jtm.command()
@click.argument('task_id')
def is_alive(task_id):
    ret = int(JtmInterface('status').call(task_id=int(task_id)))
    if ret == -88:
        eprint("jtm-isalive: command timeout.")
        sys.exit(-1)
    elif ret in [0, 1, 2]:
        print("yes")
        sys.exit(0)
    else:
        print("no")
        sys.exit(1)

@jtm.command()
@click.argument('task_id')
def kill(task_id):
    ret = int(JtmInterface('kill').call(task_id=task_id))
    if ret == -88:
        eprint("jtm-kill: command timeout.")
        sys.exit(-1)
    elif ret == -5:
        eprint("jtm-kill: task id not found.")
        sys.exit(-1)
    sys.exit(0) if ret == 0 else sys.exit(1)


@jtm.command()
@click.option('--task-pool', default=None, help="task pool")
@click.option('--cluster', type=click.Choice(COMPUTE_RESOURCES), help="cluster")
@click.option('--loglevel', default=None, help="Set loglevel (default=info)")
def remove_pool(task_pool, cluster, loglevel):
    ret = int(JtmInterface('remove_pool').call(task_pool=task_pool,
                                               jtm_host_name=cluster,
                                               log_level=loglevel))
    print(("removed" if ret == 1 else "failed"))
    sys.exit(0) if ret > 0 else sys.exit(1)


@jtm.command()
@click.argument('task_id')
def resource_log(task_id):
    resourceFile = JtmInterface('resource').call(task_id=task_id)
    if resourceFile == -88:
        eprint("jtm-resource-log: command timeout.")
        sys.exit(-1)
    if os.path.isfile(resourceFile):
        f = open(resourceFile, 'rU')
        reader = csv.DictReader(f, fieldnames=(
            "child_pid",            #1
            "clone_time_rate",      #2
            "cpu_load",             #3
            "end_date",             #4
            "host_name",            #5
            "ip_address",           #6
            "job_time",             #7
            "life_left",            #8
            "mem_per_core",         #9
            "mem_per_node",         #10
            "num_cores",            #11
            "num_tasks",            #12
            "num_workers_on_node",  #13
            "perc_mem_used",        #14
            "pool_name",            #15
            "ret_msg",              #16
            "rmem_usage",           #17
            "root_pid",             #18
            "run_time",             #19
            "slurm_jobid",          #20
            "task_id",              #21
            "vmem_usage",           #22
            "worker_id",            #23
            "worker_type",          #24
            "jtm_host_name",        #25
            "nwpn"                  #26
        ))
        print("""
    {
      "task_id": %d,
      "resource_log": %s
    }""" % (task_id, json.dumps([row for row in reader])))
    else:
        resourceFile = None
        print(("Rekource file, %s, not found." % (resourceFile)))

    sys.exit(0) if resourceFile is not None else sys.exit(1)


@jtm.command()
##todo: make this a mutually exclusive group (see https://github.com/click-contrib/click-option-group)
@click.option('--cromwell')
@click.option('--task-json-file')
@click.option('--task-json-str')
#todo end
@click.option('--cluster', type=click.Choice(COMPUTE_RESOURCES), help="cluster")
@click.option('--loglevel', default=None, help="Set loglevel (default=info)")
# resource requirement
# This will be used for create a separate pool of jtm-worker(s)
@click.option('--account', default=CORI_ACCNT)
@click.option('--cpu', type=int)
@click.option('--constraint', default=CORI_CONSTRAINT, type=click.Choice(["haswell", "knl", "skylake"]), help="Set the architecture to Haswell or KNL on Cori")
@click.option('--job-id', help="Unique Cromwell job id with step name")
@click.option('--memory', default=MEMPERNODE, help="Unique Cromwell job id with step name")
@click.option('--num-nodes', default=NNODES, type=int, help="Number of nodes for the pool. Default=1.")
@click.option('--num-workers-per-node', default=NWORKERS, type=int, help="Number of worker per node. Default=1")
@click.option('--pool-name', default='small')
@click.option('--qos', default=CORI_QOS, type=click.Choice(["genepool_special", "genepool_shared", "jgi_shared", "jgi_exvivo"]))
@click.option('--shared', default=1, type=int, help='shared workers')
@click.option('--time')
def submit(cromwell, task_json_file, task_json_str, cluster, loglevel, account, cpu, constraint, job_id, memory, nnodes, nwpn, pool, qos, shared, time):
    if time:
        assert cpu and memory and numNodes and pool, "USAGE: runtime (--time) should be set with 'cpu' (--cpu), memory (--memory), node (--num-nodes), and pool name (--pool-name)."

    if pool == "default":
        pool = "small"

    jsonTask = ""
    if cromwell:
        jsonTask = cromwell
    elif task_json_str:
        jsonTask = task_json_str

    ret = int(JtmInterface('task').call(task_file=task_json_file,
                                        task_json=jsonTask,
                                        task_id=0,
                                        jtm_host_name=cluster,
                                        job_time=time,
                                        node_mem=memory,
                                        num_core=cpu,
                                        pool_name=pool,
                                        log_level=loglevel,
                                        constraint=constraint,
                                        qos=qos,
                                        shared=shared,
                                        nwpn=nwpn,
                                        node=nnodes,
                                        job_id=job_id,
                                        account=account))

    if ret == -5:
        eprint("jtm-submit: invalid task or runtime definition.")
        return 1
    elif ret == -88:
        eprint("jtm-submit: command timeout.")
        return -1

    print(("JTM task ID %d" % ret))
    return 0 if ret != 0 else 1


@jtm.command()
@click.option('--heartbeat', type=int, default=10, help="heartbeat interval in second")
@click.option('--loglevel', default='info', help="Set loglevel (default=info)")
@click.option('--jobdir', default=JOB_LOG, help="jtm job file path")
@click.option('--logdir', default=JTM_LOG, help="jtm log file path")
@click.option('--task-pool', default='small', help="Set user-defined pool name. This should be same with task's 'pool' name.")
@click.option('--timeout', type=int, help="Set the timer for worker to terminate. If there is no request from the client for the specified seconds, the worker terminates itself (default: 60 seconds)")
@click.option('--zerofile', default=False, type=bool, help="Allow zero size output file(s).")
@click.option('--dry-run', default=False, type=bool, help="dry run")
@click.option('--job-id', default=0, type=int, help="Slurm job id")
# worker type params
@click.option('--worker-type', type=click.Choice(['static', 'dynamic']), help="jtm worker type")
@click.option('--cluster', type=click.Choice(COMPUTE_RESOURCES), help="cluster")
@click.option('--clone-time', type=float, help="cloning timing")
@click.option('--num-workers', default=1, type=int, help="number of worker per node.")
@click.option('--worker-id', help="worker id for dynamic workers")
# slurm related params
@click.option('--account', help="Charge  resources used by this job to specified account.")
@click.option('--nodes', default=NNODES, type=int, help="number of nodes.")
@click.option('--cpus-per-task', type=int, help="number of cpus")
@click.option('--constraint', default="haswell", type=click.Choice(["haswell", "knl", "skylake"]), help="Cori constraint, haswell or knl (default: haswell)")
@click.option('--mem', help="specify the real memory required per node")
@click.option('--mem-per-cpu', help="minimum memory required per allocated CPU")
@click.option('--qos', default=CORI_QOS, type=click.Choice(["genepool_special", "genepool_shared", "jgi_shared", "jgi_exvivo"]))
@click.option('--time', help="limit on the total run time")
@click.option('--job-name', help="Slurm job name")
def worker(heartbeat, loglevel, jobdir, logdir, task_pool, timeout, zerofile, dry_run, job_id, worker_type, cluster, clone_time, num_workers, worker_id, account, num_nodes, num_cpus, constraint, mem, mem_per_cpu, qos, time, job_name):
    # Note: this is for allowing zero sized output file when check if the output file is
    #  successfully created. Currently output file checking is opted out in JTM.
    global g_allowZeroOutputFile
    g_allowZeroOutputFile = zerofile

    # Set uniq worker id if worker id is provided in the params
    if worker_id:
        global g_uniqworkerid
        g_uniqworkerid = worker_id

    # Logger setting
    logger.info("Set jtm log file location to %s", logdir)
    logDir = logdir
    make_dir_p(logDir)

    # Job dir setting
    logger.info("Set jtm job file location to %s", jobdir)
    jobDir = jobdir
    make_dir_p(jobDir)

    logger.info("Unique worker ID: %s", g_uniqWorkerId)
    logger.info("**************")
    logger.info("Run mode: %s", "prod" if PRODUCTION else "dev")
    logger.info("**************")

    setup_custom_logger(loglevel, logDir, JTM_WORKER_STREAM_LOGGING, JTM_WORKER_FILE_LOGGING, workerId=g_uniqWorkerId)

    # Slurm config
    nNodes = 0
    if num_nodes:
        nNodes = num_nodes
        assert mem is not None, "-N needs --mem-per-cpu (-mc) setting."


    ###########################################################################
    # 11.13.2018 decided to remove all default values from argparse
    numWorkersPerNode = num_workers if num_workers else 1
    nCpus = num_cpus if num_cpus else 1
    memPerCpu = mem_per_cpu if mem_per_cpu else "1GB"
    memPerNode = mem if mem else ""
    ###########################################################################

    jobTime = time if time else JOBTIME
    accnt = account if account else CORI_ACCNT
    qos = qos if qos else CORI_QOS
    workerType = worker_type
    jobname = "jtm_" + job_name if job_name else "jtm_%s_worker" % (workerType)

    # Set task queue name
    innerTaskReqQ = None
    hbInterval = heartbeat if heartbeat else WORKER_HB_SEND_INTERVAL
    timeOut = timeout if timeout else WORKER_TIMEOUT  # worker timeout. Not used

    # If you want to create a custom worker pool, use this.
    # JTM detects "pool" field from task json, and creates custom pool if it's set.
    # All tasks with the pool name will be directed to the custom pool
    # userAccnt = getpass.getuser().lower()

    if task_pool:
        innerTaskReqQ = JTM_INNER_REQUEST_Q + "." + task_pool
    else:  # not reachable. default innerTaskReqQ = small
        innerTaskReqQ = JTM_INNER_REQUEST_Q

    cloneTimeRate = clone_time if clone_time else CTR
    if workerType in ("static", "dynamic"):
        assert cluster is not "" and cluster is not "local", "Static or dynamic worker needs a cluster setting (--cluster)."

    slurmJobId = job_id
    clusterName = cluster
    dryRun = dry_run

    if clusterName == "cori" and memPerCpu != "" and float(memPerCpu.replace("GB", "").replace("G", "")) > 1.0:
        logger.critical("--mem-per-cpu in Cori shouldn't be larger than 1GB. User '--mem' instead.")
        sys.exit(1)

    logger.info("RabbitMQ broker: %s", RMQ_HOST)
    logger.info("Task queue name: %s", innerTaskReqQ)
    logger.info("Worker type: %s", worker_type)

    # sbatch static worker
    # -N: If -N is not specified, the default behavior is to allocate enough nodes to satisfy
    #     the requirements of the -n and -c options. Recommended to set --mem
    # -c: Set the num cores needed. Recommended to set --mem-per-cpu

    # SBATCH
    # -N, --nodes=<minnodes[-maxnodes]>
    #               Request  that  a minimum of minnodes nodes be allocated to this job.  A maximum node count may also
    #               be specified with maxnodes.  If only one number is specified, this is used as both the minimum  and
    #               maximum  node count.  The partition's node limits supersede those of the job.  If a job's node lim-
    #               its are outside of the range permitted for its associated partition, the job  will  be  left  in  a
    #               PENDING  state.   This  permits  possible  execution  at  a later time, when the partition limit is
    #               changed.  If a job node limit exceeds the number of nodes configured in the partition, the job will
    #               be  rejected.  Note that the environment variable SLURM_JOB_NODES will be set to the count of nodes
    #               actually allocated to the job. See the ENVIRONMENT VARIABLES  section for more information.  If  -N
    #               is  not  specified, the default behavior is to allocate enough nodes to satisfy the requirements of
    #               the -n and -c options.  The job will be allocated as many nodes as possible within the range speci-
    #               fied  and  without  delaying the initiation of the job.  The node count specification may include a
    #               numeric value followed by a suffix of "k" (multiplies numeric value by 1,024)  or  "m"  (multiplies
    #               numeric value by 1,048,576).
    # -c, --cpus-per-task=<ncpus>
    #               Advise  the  Slurm  controller  that  ensuing job steps will require ncpus number of processors per
    #               task.  Without this option, the controller will just try to allocate one processor per task.
    #
    #               For instance, consider an application that has 4 tasks, each requiring 3 processors.  If our  clus-
    #               ter is comprised of quad-processors nodes and we simply ask for 12 processors, the controller might
    #               give us only 3 nodes.  However, by using the --cpus-per-task=3 options, the controller  knows  that
    #               each  task requires 3 processors on the same node, and the controller will grant an allocation of 4
    #               nodes, one for each of the 4 tasks.
    # -n, --ntasks=<number>
    #               sbatch  does  not  launch tasks, it requests an allocation of resources and submits a batch script.
    #               This option advises the Slurm controller that job steps run within the  allocation  will  launch  a
    #               maximum of number tasks and to provide for sufficient resources.  THE DEFAULT IS ONE TASK PER NODE,
    #               BUT NOTE THAT THE --cpus-per-task OPTION WILL CHANGE THIS DEFAULT.
    # --mem=<size[units]>
    #               Specify the real memory required per node.  Default units are megabytes unless the SchedulerParame-
    #               ters  configuration  parameter includes the "default_gbytes" option for gigabytes.  Different units
    #               can be specified using the suffix [K|M|G|T].  Default value is DefMemPerNode and the maximum  value
    #               is  MaxMemPerNode.  If  configured, both parameters can be seen using the scontrol show config com-
    #               mand.  THIS PARAMETER WOULD GENERALLY BE USED  IF  WHOLE  NODES  ARE  ALLOCATED  TO  JOBS  (Select-
    #               Type=select/linear).  Also see --mem-per-cpu.  --mem and --mem-per-cpu are mutually exclusive.
    #
    #               NOTE: A memory size specification of zero is treated as a special case and grants the job access to
    #               all of the memory on each node.  If the job is allocated multiple nodes in a heterogeneous cluster,
    #               the  memory  limit on each node will be that of the node in the allocation with the smallest memory
    #               size (same limit will apply to every node in the job's allocation).
    #
    #               NOTE: Enforcement of memory limits currently relies upon the task/cgroup plugin or enabling of
    #               accounting, which samples memory use on a periodic basis (data need not be stored, just collected).
    #               In both cases memory use is based upon the job's Resident Set Size (RSS). A  task  may  exceed  the
    #               memory limit until the next periodic accounting sample.
    # --mem-per-cpu=<size[units]>
    #               Minimum memory required per allocated CPU.  Default units are megabytes unless the SchedulerParame-
    #               ters configuration parameter includes the "default_gbytes" option for gigabytes.  Default value  is
    #               DefMemPerCPU  and  the  maximum  value  is  MaxMemPerCPU (see exception below). If configured, both
    #               parameters can  be  seen  using  the  scontrol  show  config  command.  Note  that  if  the  job's
    #               --mem-per-cpu value exceeds the configured MaxMemPerCPU, then the user's limit will be treated as a
    #               memory limit per task; --mem-per-cpu will be reduced  to  a  value  no  larger  than  MaxMemPerCPU;
    #               --cpus-per-task  will  be  set and the value of --cpus-per-task multiplied by the new --mem-per-cpu
    #               value will equal the original --mem-per-cpu value specified by the user.  This parameter would gen-
    #               erally  be  used  if  individual processors are allocated to jobs (SelectType=select/cons_res).  If
    #               resources are allocated by the core, socket or whole nodes; the number of CPUs allocated to  a  job
    #               may  be  higher  than the task count and the value of --mem-per-cpu should be adjusted accordingly.
    #               Also see --mem.  --mem AND --mem-per-cpu ARE MUTUALLY EXCLUSIVE.
    #
    # NOTE: Cori shared queue doesn't need "account" and "qos"
    #       If you set --qos=genepool or --qos=genepool_shared, you don't need to set "-C haswell"
    #                                                           you should set "-A fungalp"
    #       If you set -A m342, no qos needed but you need to set "-C haswell"
    #
    if slurmJobId == 0 and workerType in ["static", "dynamic"]:
        jobFile = os.path.join(jobDir, "jtm_%s_worker_%s.job" % (workerType, g_uniqWorkerId))
        jobStr = ""
        otherParams = ""

        if clusterName in ("denovo", "cori", "lawrencium"):

            with open(jobFile, "w") as jf:
                jobStr += "#!/bin/bash -l"

                if clusterName in ("cori", "denovo"):

                    if num_nodes:
                        jobStr += """
#SBATCH -N %(nnodes)d
#SBATCH --mem=%(mem)s""" % dict(nnodes=nNodes, mem=memPerNode)
                        otherParams += " -N %(nnodes)d -m %(mem)s" % dict(nnodes=nNodes, mem=memPerNode)

                        if num_cpus:
                            jobStr += """
#SBATCH -c %(ncores)d""" % dict(ncores=nCpus)
                            otherParams += " -c %(ncores)d" % dict(ncores=nCpus)

                    else:
                        jobStr += """
#SBATCH -c %(ncores)d""" % dict(ncores=nCpus)
                        otherParams += " -c %(ncores)d" % dict(ncores=nCpus)

                        if mem_per_node:
                            jobStr += """
#SBATCH --mem=%(mem)s""" % dict(mem=memPerNode)
                            otherParams += " -m %(mem)s " % dict(mem=memPerNode)
                        else:
                            jobStr += """
#SBATCH --mem-per-cpu=%(mempercore)s""" % dict(mempercore=memPerCpu)
                            otherParams += " -mc %(mempercore)s" % dict(mempercore=memPerCpu)

                        if worker_id:
                            otherParams += " -wi %(wid)s${i}" % dict(wid=g_uniqWorkerId)

                    if clusterName == "denovo":
                        pass
                    elif clusterName == "cori":

                        # Need to set both --qos=genepool (or genepool_shared) _and_ -A fungalp
                        # OR
                        # no qos _and_ -A m342 _and_ -C haswell

                        # Note: currently constraint in ["haswell" | "knl"]
                        if constraint == "haswell":
                            if qos:
                                jobStr += """
#SBATCH -q %(qosname)s""" % dict(qosname=qos)
                                otherParams += " -q %(qosname)s" % dict(qosname=qos)

                            else:
                                jobStr += """
#SBATCH -q %(qosname)s""" % dict(qosname=qos)

                            jobStr += """
#SBATCH -C haswell"""
                            if accnt == "m342":
                                otherParams += " -A %(sa)s" % dict(sa="m342")

                            jobStr += """
#SBATCH -A %(accnt)s""" % dict(accnt=accnt)

                        elif constraint == "knl":
                            # Note: Basic KNL setting = "-q regular -A m342 -C knl"
                            #
                            # Note: KNL MCDRAM setting -> cache or flat
                            #  cache mode - MCDRAM is configured entirely as a last-level cache (L3)
                            #  flat mode - MCDRAM is configured entirely as addressable memory
                            #  ex) #SBATCH -C knl,quad,cache
                            #  ex) #SBATCH -C knl,quad,flat
                            #      --> srun <srun options> numactl -p 1 yourapplication.x
                            #
                            # Note: for knl, we should use m342
                            jobStr += """
#SBATCH -C knl
#SBATCH -A m342
#SBATCH -q regular"""

                        elif constraint == "skylake":
                            jobStr += """
#SBATCH -C skylake
#SBATCH -A %(accnt)s
#SBATCH -q %(qosname)s""" % dict(accnt=accnt, qosname=qos)
                            otherParams += " -A %(accnt)s -q %(qosname)s" % dict(accnt=accnt, qosname=qos)

                        jobStr += """
#SBATCH -t %(walltime)s
#SBATCH --job-name=%(jobname)s
#SBATCH -o %(jobdir)s/jtm_%(wtype)s_worker_%(wid)s.out
#SBATCH -e %(jobdir)s/jtm_%(wtype)s_worker_%(wid)s.err
%(exclusive)s

module unload python
%(envactivation)s
for i in {1..%(numworkerspernode)d}
do
    echo "jobid: $SLURM_JOB_ID"
    jtm-worker --jobid $SLURM_JOB_ID \
-cl cori \
-wt %(wtype)s \
-t %(walltime)s \
-ct %(ctr)f %(tp)s \
-nw %(numworkerspernode)d \
-C %(constraint)s \
%(otherparams)s &
    sleep 1
done
wait
""" % dict(walltime=jobTime,
           jobdir=jobDir,
           wid=g_uniqWorkerId,
           wtype=workerType,
           ctr=cloneTimeRate,
           tp="-tp " + task_pool if task_pool else "",
           numworkerspernode=numWorkersPerNode,
           envactivation=ENV_ACTIVATION,
           otherparams=otherParams,
           constraint=constraint,
           jobname=jobname,
           exclusive="#SBATCH --exclusive" if constraint != "skylake" else "")



                elif clusterName == "lawrencium":

                    if worker_id:
                        otherParams += " -wi %s" % (worker_id)

                    jobStr += """
#SBATCH --time=%(walltime)s
#SBATCH --job-name=%(jobname)s
#SBATCH --partition=%(partitioname)s
#SBATCH --qos=%(qosname)s
#SBATCH --account=%(accnt)s
#SBATCH --nodes=%(nnodes)d
#SBATCH --mem=%(mem)s
#SBATCH -o %(jobdir)s/jtm_%(wtype)s_worker_%(wid)s.out
#SBATCH -e %(jobdir)s/jtm_%(wtype)s_worker_%(wid)s.err
%(envactivation)s

for i in {1..%(numworkerspernode)d}
do
    echo "jobid: $SLURM_JOB_ID"
    jtm-worker --jobid $SLURM_JOB_ID \
-cl lawrencium \
-wt %(wtype)s \
-t %(walltime)s \
-ct %(ctr)f %(tp)s \
-nw %(numworkerspernode)d \
%(otherparams)s &
    sleep 1
done
wait
""" % dict(walltime=jobTime,
           wtype=workerType,
           partitioname=LAWRENCIUM_PARTITION,
           qosname=LAWRENCIUM_QOS,
           accnt=LAWRENCIUM_ACCNT,
           nnodes=nNodes if num_nodes else 1,
           mem=memPerNode,
           jobdir=jobDir,
           wid=g_uniqWorkerId,
           envactivation=ENV_ACTIVATION,
           ctr=cloneTimeRate,
           tp="-tp " + task_pool if task_pool else "",
           numworkerspernode=numWorkersPerNode,
           otherparams=otherParams,
           jobname=jobname)


                jf.writelines(jobStr)

            if dryRun:
                print(jobStr)
                sys.exit(0)

            so, se, ec = run_sh_command("sbatch --parsable %s" % (jobFile), live=True, log=logger)
            assert ec == 0, "Failed to run jtm-worker to sbatch dynamic worker."
            return ec
            # if ec == 0:
            #     sys.exit(0)
            # else:
            #     logger.critical("Failed to run jtm-worker to sbatch dynamic worker.")
            #     sys.exit(1)

        elif clusterName == "aws":
            pass


    # If it's spawned by sbatch
    # TODO: need to record job_id, worker_id, worker_type, starting_time, wallclocktime
    # scontrol show jobid -dd <jobid> ==> EndTime
    # scontrol show jobid <jobid> ==> EndTime
    # sstat --format=AveCPU,AvePages,AveRSS,AveVMSize,JobID -j <jobid> --allsteps
    #
    # if endtime - starttime <= 10%, execute sbatch again
    # if slurmJobId != 0 and workerType == "static":
    #     logger.debug("workerType: {}".format(workerType))
    #     logger.debug("slurmJobId: {}".format(slurmJobId))

    # Dynamic workers creates [[two]] children when it approaches to the wallclocktime limit
    # considering the task queue length
    # Also, maintain the already requested number of workers
    # if no more workers needed, it won't call sbatch
    # elif slurmJobId != 0 and workerType == "dynamic":
    #     logger.debug("workerType: {}".format(workerType))
    #     logger.debug("slurmJobId: {}".format(slurmJobId))

    # Remote broker (mq.nersc.gov)
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()
    # ch.confirm_delivery()
    ch.exchange_declare(exchange=JTM_INNER_MAIN_EXCH,
                        exchange_type="direct",
                        passive=False,
                        durable=True,
                        auto_delete=False)

    # Declare task receiving queue (client --> worker)
    #
    # If you have a queueu that is durable, RabbitMQ will never lose our queue.
    # If you have a queue that is exclusive, then when the channel that declared
    # the queue is closed, the queue is deleted.
    # If you have a queue that is auto-deleted, then when there are no
    # subscriptions left on that queue it will be deleted.
    #
    ch.queue_declare(queue=innerTaskReqQ,
                     durable=True,
                     exclusive=False,
                     auto_delete=True)
    ch.queue_bind(exchange=JTM_INNER_MAIN_EXCH,
                  queue=innerTaskReqQ,
                  routing_key=innerTaskReqQ)

    logger.info("Waiting for a request...")

    # Start task termination thread
    global g_taskKillProc
    g_taskKillProc = multiprocessing.Process(target=recv_task_kill_request_thread,
                                             args=(task_pool,
                                                   g_uniqWorkerId,
                                                   clusterName))
    g_taskKillProc.daemon = True
    g_taskKillProc.start()

    # Start hb receive thread
    global g_recvHbFromClientProc
    g_recvHbFromClientProc = multiprocessing.Process(target=send_hb_to_client_thread,
                                                     args=(os.getpid(),
                                                           hbInterval,
                                                           g_userProcPid,  # this pid is non-zero only after user task process is forked.
                                                           g_uniqWorkerId,
                                                           taskIdRecvP,
                                                           slurmJobId,
                                                           workerType,
                                                           memPerNode,
                                                           memPerCpu,
                                                           nCpus,
                                                           jobTime,
                                                           cloneTimeRate,
                                                           innerTaskReqQ,
                                                           task_pool if task_pool else "",  # need the original task queue name
                                                           numWorkersPerNode))

    g_recvHbFromClientProc.daemon = True
    g_recvHbFromClientProc.start()
    logger.info("Start sending my heartbeat to the client in every %d sec to %s" % (hbInterval, WORKER_HB_Q_POSTFIX))

    # Start poison receive thread
    global g_recvPoisonProc
    g_recvPoisonProc = multiprocessing.Process(target=recv_reproduce_or_die_thread,
                                               args=(task_pool,  # need the original task queue name
                                                     g_uniqWorkerId,
                                                     clusterName,
                                                     memPerNode,
                                                     memPerCpu,
                                                     nNodes,
                                                     nCpus,
                                                     jobTime,
                                                     cloneTimeRate,
                                                     numWorkersPerNode))
    g_recvPoisonProc.daemon = True
    g_recvPoisonProc.start()

    # Start hb send thread
    global g_sendHbToClientProc
    g_sendHbToClientProc = multiprocessing.Process(target=recv_hb_from_client_thread,
                                                   args=(innerTaskReqQ,
                                                         timeOut,
                                                         os.getpid(),
                                                         g_uniqWorkerId))
    g_sendHbToClientProc.daemon = True
    g_sendHbToClientProc.start()

    if timeOut != 0:
        logger.info("The worker timeout is set to %s sec. Will not be terminated even without jtm's heartbeat.", timeOut)

    # Waiting for request
    ch.basic_qos(prefetch_count=1)

    # OLD
    if int(PIKA_VER[0]) < 1:  # v0.13.1
        ch.basic_consume(on_request, queue=innerTaskReqQ, no_ack=False)
    else:  # v1.0.1 or higher
        ch.basic_consume(queue=innerTaskReqQ, on_message_callback=on_request, auto_ack=False)

    # NEW
    # Ref) https://github.com/pika/pika/blob/1.0.1/examples/basic_consumer_threaded.py
    #      https://stackoverflow.com/questions/51752890/how-to-disable-heartbeats-with-pika-and-rabbitmq
    # threads = []
    # on_message_callback = functools.partial(on_task_request, args=(conn, threads))
    # if int(PIKA_VER[0]) < 1:  # v0.13.1
    #     ch.basic_consume(innerTaskReqQ, on_message_callback)
    # else:  # v1.0.1 or higher
    #     ch.basic_consume(queue=innerTaskReqQ, on_message_callback=on_message_callback)

    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        if g_recvHbFromClientProc: g_recvHbFromClientProc.terminate()
        if g_sendHbToClientProc: g_sendHbToClientProc.terminate()
        if g_recvPoisonProc: g_recvPoisonProc.terminate()
        if g_taskKillProc: g_taskKillProc.terminate()
        ch.stop_consuming()

    # Wait for all to complete
    # Note: prefetch_count=1 ==> #thread = 1
    # for thread in threads:
    #     thread.join()

    # Unreachable
    if ch: ch.close()
    if conn: conn.close()


    return 0


@jtm.command()
@click.option('--loglevel', default='info', help="Set loglevel")
@click.option('--logdir', default=JTM_LOG, help="jtm log file path")
@click.option('--resource', default=False, type=bool, help="Print resource usage report. Display format: ['child_pid', 'clone_time_rate', 'cpu_load', 'end_date', 'host_name', 'ip_address', 'job_time', 'life_left', 'mem_per_core', 'mem_per_node', 'num_cores', 'num_tasks', 'num_workers_on_this_node', 'perc_mem_used', 'pool_name', 'ret_msg', 'rmem_usage', 'root_pid', 'run_time', 'slurm_jobid', 'task_id', 'vmem_usage', 'worker_id', 'worker_type'].")
def manager(loglevel, logdir, resource):
    # Logger setting
    logger.info("Set jtm log file location to %s", logdir)
    logDir = logdir
    make_dir_p(logDir)

    setup_custom_logger(loglevel, logdir, JTM_STREAM_LOGGING, JTM_FILE_LOGGING)

    # Remote broker (rmq.nersc.gov) connection open
    rabbitConnection = RmqConnectionHB()
    conn = rabbitConnection.open()
    ch = conn.channel()
    # ch.confirm_delivery()  # 11192018 to test task is not discarded
    ch.exchange_declare(exchange=JGI_JTM_MAIN_EXCH,
                        exchange_type="direct",
                        durable=True,
                        auto_delete=False)
    logger.info("JTM main exchange: %s", JGI_JTM_MAIN_EXCH)

    # Declare task sending queue (client --> worker)
    #
    # This queue can be declared from worker side. BUT
    # This should be done here to enable client to send the tasks to the
    # task queue. If it is not called here, the client should be run after
    # checking if the worker is already running so that we can ensure that
    # the task queue is already declared and safe to be used by the client.
    #
    try:
        # exclusive=False -> do not remove the queue even when the connection is closed.
        ch.queue_declare(queue=JTM_TASK_REQUEST_Q,
                         durable=True,
                         exclusive=False,
                         auto_delete=True)
        ch.queue_declare(queue=JTM_TASK_RESULT_Q,
                         durable=True,
                         exclusive=False,
                         auto_delete=True)
    except Exception as detail:
        logger.exception("Exception: The queue, %s is already in use.", JTM_TASK_RESULT_Q)
        logger.exception("Detail: %s", str(detail))
        # sys.exit(1)
        return 1

    # Queue binding for getting task request from JAWS
    ch.queue_bind(exchange=JGI_JTM_MAIN_EXCH,
                  queue=JTM_TASK_REQUEST_Q,
                  routing_key=JTM_TASK_REQUEST_Q)

    logger.info("RabbitMQ broker: %s", RMQ_HOST)
    logger.info("RabbitMQ port: %s", RMQ_PORT)
    logger.info("Default task queue name: %s", JTM_TASK_REQUEST_Q)
    logger.info("Default result queue name: %s", JTM_TASK_RESULT_Q)
    logger.info("Pika version: %s", PIKA_VER)
    logger.info("Database server: %s", MYSQL_HOST)
    logger.info("Database name: %s", MYSQL_DB)
    logger.info("JTM user name: %s", USER_NAME)
    logger.info("**************")
    logger.info("Run mode: %s", "prod" if PRODUCTION else "dev")
    logger.info("**************")

    #
    # MySQL: prepare task table
    #
    db = DbSqlMy(db=MYSQL_DB)
    db.ddl(JTM_SQL["set_timezone"])
    db.ddl(JTM_SQL["create_database"] % MYSQL_DB)
    db.ddl(JTM_SQL["use_database"] % MYSQL_DB)
    db.ddl(JTM_SQL["create_table_tasks"])
    db.ddl(JTM_SQL["create_table_runs"])
    # db.createIndices(names=["taskId"], table="runs")
    db.ddl(JTM_SQL["create_table_workers"])
    db.close()


    # Start heartbeat checking thread
    sendHbToWorkersProc = multiprocessing.Process(target=send_hb_to_workers_thread)
    sendHbToWorkersProc.daemon = True
    sendHbToWorkersProc.start()
    logger.info("Broadcasting heartbeats to workers...")

    # Start heartbeat receiving thread
    wkHbQueueName = WORKER_HB_Q_POSTFIX
    recvHbFromWorkersProc = multiprocessing.Process(target=recv_hb_from_workers_thread,
                                                    args=(wkHbQueueName,
                                                          os.getpid(),
                                                          sendHbToWorkersProc,
                                                          logdir,
                                                          resource))
    recvHbFromWorkersProc.daemon = True
    recvHbFromWorkersProc.start()
    logger.info("Waiting for worker\'s heartbeats from %s", wkHbQueueName)

    # Start task kill thread
    taskKillProc = multiprocessing.Process(target=task_kill_thread)
    taskKillProc.daemon = True
    taskKillProc.start()

    # Start worker cleanup thread
    workerCleanupProc = multiprocessing.Process(target=zombie_worker_cleanup_thead)
    workerCleanupProc.daemon = True
    workerCleanupProc.start()

    # Start result receiving thread
    # listening to JTM_INNER_RESULT_Q to which all workers will send result messages
    recvResFromWorkersProc = multiprocessing.Process(target=recv_result_from_workers_thread)
    recvResFromWorkersProc.daemon = True
    recvResFromWorkersProc.start()
    logger.info("Waiting for a task request from %s", JTM_TASK_REQUEST_Q)

    # basic_qos():
    # prefetch_size (int) – This field specifies the prefetch window size. The server will send a message in advance
    # if it is equal to or smaller in size than the available prefetch size (and also falls into other prefetch limits).
    # May be set to zero, meaning “no specific limit”, although other prefetch limits may still apply. The prefetch-size
    # is ignored if the no-ack option is set in the consumer.
    # prefetch_count (int) – Specifies a prefetch window in terms of whole messages. This field may be used in
    # combination with the prefetch-size field; a message will only be sent in advance if both prefetch windows (and
    # those at the channel and connection level) allow it. The prefetch-count is ignored if the no-ack option is set in
    # the consumer.
    # all_channels (bool) – Should the QoS apply to all channels
    #
    ch.basic_qos(prefetch_count=1)
    if int(PIKA_VER[0]) < 1:  # v0.13.1
        ch.basic_consume(on_task_request, queue=JTM_TASK_REQUEST_Q, no_ack=False)
    else:  # v1.0.1 or higher
        ch.basic_consume(queue=JTM_TASK_REQUEST_Q, on_message_callback=on_task_request, auto_ack=False)

    # NOTE: the below methods might cause error in consuming messages which are already queued
    # ch.basic_consume(lambda ch, method, properties, body: on_task_request(ch, method, properties, body, g_dbConnPool),
    #                  queue=JTM_TASK_REQUEST_Q,
    #                  auto_ack=False)
    # on_task_request_callback = functools.partial(on_task_request, args=(g_dbConnPool))
    # ch.basic_consume(on_task_request_callback,
    #                  JTM_TASK_REQUEST_Q,
    #                  auto_ack=False)

    # Keep consuming messages from task request queue
    try:
        ch.start_consuming()
    except KeyboardInterrupt:
        if sendHbToWorkersProc: sendHbToWorkersProc.terminate()
        if recvHbFromWorkersProc: recvHbFromWorkersProc.terminate()
        if recvResFromWorkersProc: recvResFromWorkersProc.terminate()
        if taskKillProc: taskKillProc.terminate()
        if workerCleanupProc: workerCleanupProc.terminate()
        ch.stop_consuming()

    # unreachable
    if ch: ch.close()
    if conn: conn.close()


    return 0
