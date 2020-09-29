#!/bin/bash

cat << EOM > /jtm/jaws_jtm/jaws-jtm.conf
[SITE]
jtm_host_name = $JAWS_SITE
user_name = jtm_$DEPLOYMENT_NAME
instance_name = \${jtm_host_name}.\${user_name}
scratch = $SITE_JTM_SCRATCH_DIR
debug = 1
[RMQ]
user = jaws
password = $SITE_RMQ_PW
host = $SITE_RMQ_HOST
port = $SITE_RMQ_PORT
vhost = jtm_$DEPLOYMENT_NAME
[SITE_RPC_CLIENT]
user = jaws
password = $SITE_RMQ_PW
host = $SITE_RMQ_HOST
port = $SITE_RMQ_PORT
vhost = jaws_$DEPLOYMENT_NAME
queue = SITE_$JAWS_SITE
[MYSQL]
host = $SITE_DB_HOST
port = $SITE_DB_PORT
user = jaws
password = $SITE_DB_PW
db = jtm_${JAWS_SITE,,}_$DEPLOYMENT_NAME
[SLURM]
# number of nodes for the pool
nnodes = 1
# num cores for small task
ncpus = 1
# default wall clock time
jobtime = 00:30:00
mempercpu = 1gb
mempernode = 5gb
partition = $SITE_CLUSTER_PARTITION
qos = $SITE_CLUSTER_QOS
charge_accnt = $SITE_CLUSTER_ACCOUNT
constraint = $SITE_CLUSTER_CONSTRAINT
[JTM]
cluster = \${SITE:jtm_host_name}
run_mode = dev
# if STANDALONE == 1, do not sending task status to JAWS Site
standalone = 0
# static worker self cloning time rate. not used.s
clone_time_rate = 0.2
# number of workers per node
num_workers_per_node = 1
log_dir = \${SITE:scratch}/jtm
# config file location of worker
worker_config_file = $SITE_JTM_WORKER_INSTALL_DIR/jaws-jtm.conf
env_activation = source $SITE_JTM_WORKER_INSTALL_DIR/bin/activate
pool_name = small
# Worker setting
# worker's hb expiration (ms)
worker_s_hb_expiration = 60000
# worker->client heartbeat sending interval in worker
worker_hb_send_interval = 2.0
# timeout for waiting the client's heartbeat; 0=no timeout. OBSOLETE
worker_timeout = 0
# zombie worker checking interval
worker_kill_interval = 300.0
# number of procs checking interval
num_procs_check_interval = 60.0
# manager's interval to check if there is task to terminate
task_kill_interval = 5.0
# client->worker heartbeat receiving interval in worker.
worker_hb_recv_interval = 5.0
# tasks status update interval. OBSOLETE
task_stat_update_interval = 0.5
# to ensure updating runs table
runs_info_update_wait = 0.5
# to ensure live worker info is updated in workers' table
worker_info_update_wait = 0.5
# interval to check the result from jtm
result_receive_interval = 0.1
# num threads to recv results from workers
num_result_recv_threads = 1
# client->worker heartbeat sending interval
client_hb_send_interval = 3.0
# live worker checking count limit
worker_hb_check_max_count = 0
# worker->client heartbeat receiving interval
client_hb_recv_interval = 6.0
# max number of trial for checking
file_checking_max_trial = 3
# sleep time between output file checking before retrial
file_check_interval = 3.0
# increase amount of wait time for file checking
file_check_int_inc = 1.5
# jtm_* CLI command timeout in secs.
jtminterface_max_trial = 3600
# worker lifeleft for task timeout
task_kill_timeout_minute = 3
jtm_inner_request_q = _jtm_inner_request_queue.\${SITE:instance_name}
jtm_inner_result_q = _jtm_inner_result_queue.\${SITE:instance_name}
jtm_task_kill_q = _jtm_task_kill.\${SITE:instance_name}
jtm_worker_poison_q = _jtm_worker_poison.\${SITE:instance_name}
jtm_task_result_q = _jtm_task_result_queue.\${SITE:instance_name}
jtm_task_request_q = _jtm_task_request_queue.\${SITE:instance_name}
jtm_status_result_q = _jtm_status_result_queue.\${SITE:instance_name}
jtm_status_request_q = _jtm_status_request_queue.\${SITE:instance_name}
worker_hb_q_postfix = _jtm_worker_hb_queue.\${SITE:instance_name}
client_hb_q_postfix = _jtm_client_hb_queue.\${SITE:instance_name}
jgi_jtm_main_exch = jgi_jtm_main_exch_\${JTM:run_mode}_\${SITE:instance_name}
jtm_client_hb_exch = jgi_jtm_client_hb_exch_\${JTM:run_mode}_\${SITE:instance_name}
jtm_inner_main_exch = jgi_jtm_inner_main_exch_\${JTM:run_mode}_\${SITE:instance_name}
jtm_worker_hb_exch = jgi_jtm_worker_hb_exch_\${JTM:run_mode}_\${SITE:instance_name}
jtm_task_kill_exch = jgi_jtm_task_kill_exch_\${JTM:run_mode}_\${SITE:instance_name}
jtm_worker_poison_exch = jgi_jtm_poison_exch_\${JTM:run_mode}_\${SITE:instance_name}
EOM


cat << EOM > /jtm/jaws_jtm_worker/jaws-jtm.conf
[SITE]
jtm_host_name = $JAWS_SITE
user_name = jtm_$DEPLOYMENT_NAME
instance_name = \${jtm_host_name}.\${user_name}
scratch = $SITE_JTM_SCRATCH_DIR
debug = 1
[RMQ]
user = jaws
password = $SITE_RMQ_PW
host = $SITE_RMQ_HOST
port = $SITE_RMQ_PORT
vhost = jtm_$DEPLOYMENT_NAME
[SITE_RPC_CLIENT]
user = jaws
password = $SITE_RMQ_PW
host = $SITE_RMQ_HOST
port = $SITE_RMQ_PORT
vhost = jaws_$DEPLOYMENT_NAME
queue = SITE_$JAWS_SITE
[MYSQL]
host = $SITE_DB_HOST
port = $SITE_DB_PORT
user = jaws
password = $SITE_DB_PW
db = jtm_${JAWS_SITE,,}_$DEPLOYMENT_NAME
[SLURM]
# number of nodes for the pool
nnodes = 1
# num cores for small task
ncpus = 1
# default wall clock time
jobtime = 00:30:00
mempercpu = 1gb
mempernode = 5gb
partition = $SITE_CLUSTER_PARTITION
qos = $SITE_CLUSTER_QOS
charge_accnt = $SITE_CLUSTER_ACCOUNT
constraint = $SITE_CLUSTER_CONSTRAINT
[JTM]
cluster = \${SITE:jtm_host_name}
run_mode = dev
# if STANDALONE == 1, do not sending task status to JAWS Site
standalone = 0
# static worker self cloning time rate. not used.s
clone_time_rate = 0.2
# number of workers per node
num_workers_per_node = 1
log_dir = \${SITE:scratch}/jtm
# config file location of worker
worker_config_file = $SITE_JTM_WORKER_INSTALL_DIR/jaws-jtm.conf
env_activation = source $SITE_JTM_WORKER_INSTALL_DIR/bin/activate
pool_name = small
# Worker setting
# worker's hb expiration (ms)
worker_s_hb_expiration = 60000
# worker->client heartbeat sending interval in worker
worker_hb_send_interval = 2.0
# timeout for waiting the client's heartbeat; 0=no timeout. OBSOLETE
worker_timeout = 0
# zombie worker checking interval
worker_kill_interval = 300.0
# number of procs checking interval
num_procs_check_interval = 60.0
# manager's interval to check if there is task to terminate
task_kill_interval = 5.0
# client->worker heartbeat receiving interval in worker.
worker_hb_recv_interval = 5.0
# tasks status update interval. OBSOLETE
task_stat_update_interval = 0.5
# to ensure updating runs table
runs_info_update_wait = 0.5
# to ensure live worker info is updated in workers' table
worker_info_update_wait = 0.5
# interval to check the result from jtm
result_receive_interval = 0.1
# num threads to recv results from workers
num_result_recv_threads = 1
# client->worker heartbeat sending interval
client_hb_send_interval = 3.0
# live worker checking count limit
worker_hb_check_max_count = 0
# worker->client heartbeat receiving interval
client_hb_recv_interval = 6.0
# max number of trial for checking
file_checking_max_trial = 3
# sleep time between output file checking before retrial
file_check_interval = 3.0
# increase amount of wait time for file checking
file_check_int_inc = 1.5
# jtm_* CLI command timeout in secs.
jtminterface_max_trial = 3600
# worker lifeleft for task timeout
task_kill_timeout_minute = 3
jtm_inner_request_q = _jtm_inner_request_queue.\${SITE:instance_name}
jtm_inner_result_q = _jtm_inner_result_queue.\${SITE:instance_name}
jtm_task_kill_q = _jtm_task_kill.\${SITE:instance_name}
jtm_worker_poison_q = _jtm_worker_poison.\${SITE:instance_name}
jtm_task_result_q = _jtm_task_result_queue.\${SITE:instance_name}
jtm_task_request_q = _jtm_task_request_queue.\${SITE:instance_name}
jtm_status_result_q = _jtm_status_result_queue.\${SITE:instance_name}
jtm_status_request_q = _jtm_status_request_queue.\${SITE:instance_name}
worker_hb_q_postfix = _jtm_worker_hb_queue.\${SITE:instance_name}
client_hb_q_postfix = _jtm_client_hb_queue.\${SITE:instance_name}
jgi_jtm_main_exch = jgi_jtm_main_exch_\${JTM:run_mode}_\${SITE:instance_name}
jtm_client_hb_exch = jgi_jtm_client_hb_exch_\${JTM:run_mode}_\${SITE:instance_name}
jtm_inner_main_exch = jgi_jtm_inner_main_exch_\${JTM:run_mode}_\${SITE:instance_name}
jtm_worker_hb_exch = jgi_jtm_worker_hb_exch_\${JTM:run_mode}_\${SITE:instance_name}
jtm_task_kill_exch = jgi_jtm_task_kill_exch_\${JTM:run_mode}_\${SITE:instance_name}
jtm_worker_poison_exch = jgi_jtm_poison_exch_\${JTM:run_mode}_\${SITE:instance_name}
EOM


# This was taken from https://github.com/giovtorres/docker-centos7-slurm/blob/master/docker-entrypoint.sh
# to avoid overwriting the original docker-entrypoint.sh which had logic to start Slurm. 
if [ ! -f "/var/lib/mysql/ibdata1" ]; then
    echo "- Initializing database"
    /usr/bin/mysql_install_db --force &> /dev/null

    echo "- Updating MySQL directory permissions"
    chown -R mysql:mysql /var/lib/mysql
    chown -R mysql:mysql /var/run/mariadb
fi

if [ ! -d "/var/lib/mysql/slurm_acct_db" ]; then
    /usr/bin/mysqld_safe --datadir="/var/lib/mysql" &

    for count in {30..0}; do
        if echo "SELECT 1" | mysql &> /dev/null; then
            break
        fi
        echo "- Starting MariaDB to create Slurm account database"
        sleep 1
    done

    if [[ "$count" -eq 0 ]]; then
        echo >&2 "MariaDB did not start"
        exit 1
    fi

    echo "- Creating Slurm acct database"
    mysql -NBe "CREATE DATABASE slurm_acct_db"
    mysql -NBe "CREATE USER 'slurm'@'localhost'"
    mysql -NBe "SET PASSWORD for 'slurm'@'localhost' = password('password')"
    mysql -NBe "GRANT USAGE ON *.* to 'slurm'@'localhost'"
    mysql -NBe "GRANT ALL PRIVILEGES on slurm_acct_db.* to 'slurm'@'localhost'"
    mysql -NBe "FLUSH PRIVILEGES"
    echo "- Slurm acct database created. Stopping MariaDB"
    killall mysqld
    for count in {30..0}; do
        if echo "SELECT 1" | mysql &> /dev/null; then
            sleep 1
        else
            break
        fi
    done
    if [[ "$count" -eq 0 ]]; then
        echo >&2 "MariaDB did not stop"
        exit 1
    fi
fi

chown slurm:slurm /var/spool/slurmd /var/run/slurmd /var/lib/slurmd /var/log/slurm

echo "- Starting all Slurm processes under supervisord"
/usr/bin/supervisord --configuration /etc/supervisord.conf

source /jtm/jaws_jtm/bin/activate && $1
