export JAWS_INSTALL_BASEDIR="/global/home/groups-sw/lr_jgicloud/jaws-install"
export JAWS_GLOBUS_EP="e8b18c38-36cd-11eb-b54c-02d9497ca481"
export JAWS_GLOBUS_HOST_PATH="/global/scratch/users/jaws"
export JAWS_LOAD_PYTHON="module load python/3.9.12"
export JAWS_PYTHON="python"
export JAWS_GROUP="jaws"
export JAWS_USERS_GROUP="jgi"
export JAWS_SCRATCH_BASEDIR="/global/scratch/users/jaws"
export JAWS_REF_DATA_DIR="/global/scratch/users/jaws/refdata"
export JAWS_MAX_RAM_GB=492
export JAWS_MAX_CPU=32
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
export JAWS_QUEUE_WAIT_SMALL="sbatch --test-only --time=72:00:00 --mem=46G --partition=lr3 --qos=condo_jgicloud --account=lr_jgicloud --wrap 'sleep 30' 2>&1"
export JAWS_QUEUE_WAIT_MEDIUM="sbatch --test-only --time=72:00:00 --mem=236G --partition=jgi --qos=normal --account=jgi --cpus-per-task=32 --wrap 'sleep 30' 2>&1"
export JAWS_QUEUE_WAIT_LARGE="sbatch --test-only --time=72:00:00 --mem=492G --partition=lr3 --qos=condo_jgicloud --constraint=lr3_c32 --account=lr_jgicloud  --wrap 'sleep 30' 2>&1"
