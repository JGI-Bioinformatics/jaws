from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.launchers import SrunLauncher, SingleNodeLauncher
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_interface, address_by_hostname
from parsl.monitoring.monitoring import MonitoringHub

cori_hsw_opts = '#SBATCH -C haswell'
cori_hsw_opts += '\n#SBATCH -A m342'

cori_gp_opts = '#SBATCH -A genother'

cori_ev_opts = '#SBATCH -M escori'
cori_ev_opts += '\n#SBATCH -A genother'

lbl_sched_opts = '#SBATCH -A jgicloud'

tahoma_sched_opts = ''

config_cori = Config(
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
                scheduler_options=cori_hsw_opts,
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
                scheduler_options=cori_gp_opts,
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
            address=address_by_interface('bond0.144'),
            provider=SlurmProvider(
                'jgi_exvivo',  # Partition / QOS
                nodes_per_block=1,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=cori_ev_opts,
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

config_lbl = Config(
    executors=[
        HighThroughputExecutor(
            label='lbl',
            address=address_by_hostname(),
            provider=SlurmProvider(
                partition='jgi',  # Partition / QOS
                nodes_per_block=1,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=lbl_sched_opts,
                # *** NOTE *** worker_init must be setup
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                # worker_init='module load python; source activate parsl-env',
                launcher=SingleNodeLauncher(),
                walltime='48:00:00',
                cmd_timeout=120,
                exclusive=False
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

config_tahoma = Config(
    executors=[
        HighThroughputExecutor(
            label='tahoma',
            address=address_by_hostname(),
            provider=SlurmProvider(
                partition='normal',  # Partition / QOS
                # https://www.emsl.pnnl.gov/MSC/UserGuide/tahoma/compute_jobs.html
                nodes_per_block=1,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=tahoma_sched_opts,
                # *** NOTE *** worker_init must be setup
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                # worker_init='module load python; source activate parsl-env',
                launcher=SingleNodeLauncher(),
                walltime='48:00:00',
                cmd_timeout=120,
                exclusive=False
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
