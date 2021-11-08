from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.launchers import SrunLauncher
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_interface, address_by_hostname
from parsl.monitoring.monitoring import MonitoringHub


CORI_HSW_OPTS = '#SBATCH -C haswell'
CORI_HSW_OPTS += '\n#SBATCH -A m342'
CORI_GP_OPTS = '#SBATCH -A genother'
CORI_EV_OPTS = '#SBATCH -M escori'
CORI_EV_OPTS += '\n#SBATCH -A genother'


CONFIG_CORI = Config(
    executors=[
        HighThroughputExecutor(
            label='cori_haswell',
            # This is the network interface on the login node to
            # which compute nodes can communicate
            address=address_by_interface('bond0.144'),
            provider=SlurmProvider(
                'regular',  # Partition / QOS
                nodes_per_block=1,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=CORI_HSW_OPTS,
                # *** NOTE *** worker_init must be setup
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                # worker_init='module load python; source activate parsl-env',
                # We request all hyperthreads on a node.
                launcher=SrunLauncher(),
                walltime='48:00:00',
                # Slurm scheduler on Cori can be slow at times,
                # increase the command timeouts
                cmd_timeout=120,
            ),
        ),
        HighThroughputExecutor(
            label='cori_genepool',
            # This is the network interface on the login node to
            # which compute nodes can communicate
            address=address_by_interface('bond0.144'),
            provider=SlurmProvider(
                'genepool',  # Partition / QOS
                nodes_per_block=1,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=CORI_GP_OPTS,
                # *** NOTE *** worker_init must be setup
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                # worker_init='module load python; source activate parsl-env',
                # We request all hyperthreads on a node.
                launcher=SrunLauncher(),
                walltime='48:00:00',
                # Slurm scheduler on Cori can be slow at times,
                # increase the command timeouts
                cmd_timeout=120,
            ),
        ),
        HighThroughputExecutor(
            label='cori_exvivo',
            # This is the network interface on the login node to
            # which compute nodes can communicate
            address=address_by_hostname(),
            provider=SlurmProvider(
                'jgi_exvivo',  # Partition / QOS
                nodes_per_block=1,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=CORI_EV_OPTS,
                # *** NOTE *** worker_init must be setup
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                # worker_init='module load python; source activate parsl-env',
                # We request all hyperthreads on a node.
                launcher=SrunLauncher(),
                walltime='48:00:00',
                # Slurm scheduler on Cori can be slow at times,
                # increase the command timeouts
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
