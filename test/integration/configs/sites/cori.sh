export JAWS_INSTALL_BASEDIR="/global/cfs/projectdirs/jaws/jaws-install"
export JAWS_GLOBUS_EP="be1ff650-dcbc-11ea-85a2-0e1702b77d41"
export JAWS_GLOBUS_HOST_PATH="/"
export JAWS_LOAD_PYTHON="module load python/3.9-anaconda-2021.11"
export JAWS_PYTHON="python3"
export JAWS_GROUP="jaws"
export JAWS_USERS_GROUP="genome"
export JAWS_SCRATCH_BASEDIR="/global/cscratch1/sd/jaws"
export JAWS_REF_DATA_DIR="/global/dna/shared/databases/jaws/refdata"
export JAWS_MAX_RAM_GB=1450
export JAWS_MAX_CPU=72
export JAWS_SUPERVISOR_PORT_PROD=64103
export JAWS_AUTH_PORT_PROD=3003
export JAWS_REST_PORT_PROD=5003
export JAWS_CROMWELL_PORT_PROD=50103
export JAWS_SUPERVISOR_PORT_STAGING=64102
export JAWS_AUTH_PORT_STAGING=3002
export JAWS_REST_PORT_STAGING=5002
export JAWS_CROMWELL_PORT_STAGING=50102
export JAWS_SUPERVISOR_PORT_DEV=64101
export JAWS_AUTH_PORT_DEV=3001
export JAWS_REST_PORT_DEV=5001
export JAWS_CROMWELL_PORT_DEV=50101
export JAWS_QUEUE_WAIT_MEDIUM="sbatch --test-only -q genepool_special -C haswell -A fungalp --time=72:00:00 --wrap 'sleep 30' 2>&1"
export JAWS_QUEUE_WAIT_XLARGE="module load esslurm && sbatch --test-only -q jgi_exvivo -C skylake -A fungalp --time=72:00:00 --wrap 'sleep 30' 2>&1"
export JAWS_FILE_SYNC_DELAY_SEC=180
