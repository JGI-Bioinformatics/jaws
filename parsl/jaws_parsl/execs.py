from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.launchers import SrunLauncher
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_interface

cori_sched_opts = '#SBATCH -C haswell'
cori_sched_opts += '\n#SBATCH -A m342'

lbl_sched_opts = '#SBATCH -A jgicloud'

cascade_sched_opts = ''

config = Config(
    executors=[
        HighThroughputExecutor(
            label='CORI',
            # This is the network interface on the login node to
            # which compute nodes can communicate
            address=address_by_interface('bond0.144'),
            provider=SlurmProvider(
                # 'regular',  # Partition / QOS
                'debug',  # Partition / QOS
                nodes_per_block=1,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=cori_sched_opts,
                # *** NOTE *** worker_init must be setup
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                worker_init='module load python; source activate parsl-env',
                # We request all hyperthreads on a node.
                launcher=SrunLauncher(),
                walltime='48:00:00',
                # Slurm scheduler on Cori can be slow at times,
                # increase the command timeouts
                cmd_timeout=120,
            ),
        ),
        HighThroughputExecutor(
            label='JGI',
            # This is the network interface on the login node to
            # which compute nodes can communicate
            address=address_by_interface('bond0.144'),
            provider=SlurmProvider(
                # 'regular',  # Partition / QOS
                'jgi',  # Partition / QOS
                nodes_per_block=2,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=lbl_sched_opts,
                # *** NOTE *** worker_init must be setup
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                worker_init='module load python; source activate parsl-env',
                # We request all hyperthreads on a node.
                launcher=SrunLauncher(),
                walltime='48:00:00',
                # Slurm scheduler on Cori can be slow at times,
                # increase the command timeouts
                cmd_timeout=120,
            ),
        ),
        HighThroughputExecutor(
            label='CASCADE',
            # This is the network interface on the login node to
            # which compute nodes can communicate
            address=address_by_interface('bond0.144'),
            provider=SlurmProvider(
                # 'regular',  # Partition / QOS
                'debug',  # Partition / QOS
                nodes_per_block=2,
                init_blocks=1,
                # string to prepend to #SBATCH blocks in the submit
                # script to the scheduler eg: '#SBATCH --constraint=knl,quad,cache'
                scheduler_options=cascade_sched_opts,
                # *** NOTE *** worker_init must be setup
                # Command to be run before starting a worker, such as:
                # 'module load Anaconda; source activate parsl_env'.
                worker_init='module load python; source activate parsl-env',
                # We request all hyperthreads on a node.
                launcher=SrunLauncher(),
                walltime='48:00:00',
                # Slurm scheduler on Cori can be slow at times,
                # increase the command timeouts
                cmd_timeout=120,
            ),
        ),
    ]
)
