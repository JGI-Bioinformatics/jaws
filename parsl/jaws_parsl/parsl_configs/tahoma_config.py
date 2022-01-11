from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.launchers import SingleNodeLauncher
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_hostname
from parsl.monitoring.monitoring import MonitoringHub

TAHOMA_SCHED_OPTS = '#SBATCH -A mscjgi'


CONFIG_TAHOMA = Config(
    executors=[
        HighThroughputExecutor(
            label='tahoma_normal',
            address=address_by_hostname(),
            provider=SlurmProvider(
                partition='normal',  # Partition / QOS
                # https://www.emsl.pnnl.gov/MSC/UserGuide/tahoma/compute_jobs.html
                nodes_per_block=1,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=TAHOMA_SCHED_OPTS,
                # *** NOTE *** worker_init must be setup
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                # worker_init='module load python; source activate parsl-env',
                launcher=SingleNodeLauncher(),
                walltime='48:00:00',
                cmd_timeout=120,
            ),
        ),
        HighThroughputExecutor(
            label='tahoma_analysis',
            address=address_by_hostname(),
            provider=SlurmProvider(
                partition='analysis',  # Partition / QOS
                # https://www.emsl.pnnl.gov/MSC/UserGuide/tahoma/compute_jobs.html
                nodes_per_block=1,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=TAHOMA_SCHED_OPTS,
                # *** NOTE *** worker_init must be setup
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                # worker_init='module load python; source activate parsl-env',
                launcher=SingleNodeLauncher(),
                walltime='48:00:00',
                cmd_timeout=120,
            ),
        ),
    ],
    monitoring=MonitoringHub(
        hub_address=address_by_hostname(),
        hub_port=55055,
        monitoring_debug=False,
        resource_monitoring_interval=10,
    )
)
