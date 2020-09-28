#!/usr/bin/env bash

cat << EOM > /etc/jaws-site.conf
[DB]
dialect = mysql+mysqlconnector
host = $DB_HOST
port = $DB_PORT
user = jaws
password = $DB_PW
db = jaws_${JAWS_SITE,,}_$DEPLOYMENT_NAME
[LOCAL_RPC_SERVER]
user = jaws
password = $RMQ_PW
host = $RMQ_HOST
port = $RMQ_PORT
vhost = jaws_$DEPLOYMENT_NAME
queue = SITE_$JAWS_SITE
num_threads = 5
max_retries = 3
[CENTRAL_RPC_SERVER]
user = jaws
password = $JAWS_CENTRAL_RMQ_PW
host = $JAWS_CENTRAL_RMQ_HOST
port = $JAWS_CENTRAL_RMQ_PORT
vhost = jaws_$DEPLOYMENT_NAME
queue = JAWS_$JAWS_SITE
num_threads = 5
max_retries = 3
[CENTRAL_RPC_CLIENT]
user = jaws
password = $JAWS_CENTRAL_RMQ_PW
host = $JAWS_CENTRAL_RMQ_HOST
port = $JAWS_CENTRAL_RMQ_PORT
vhost = jaws_$DEPLOYMENT_NAME
queue = CENTRAL
[GLOBUS]
client_id = $JAWS_GLOBUS_CLIENT_ID
endpoint_id = $GLOBUS_EP
root_dir = $GLOBUS_ROOT_DIR
default_dir = $GLOBUS_DEFAULT_DIR
[SITE]
id = $JAWS_SITE
uploads_subdirectory = $UPLOADS_SUBDIR
downloads_subdirectory = $DOWNLOADS_SUBDIR
[CROMWELL]
url = http://localhost:$CROMWELL_PORT
EOM

SUPERVISOR_JAWS_CONF=""
cat << EOM > /etc/supervisord.conf
[supervisord]
logfile=/tmp/jaws-log
pidfile=/tmp/jaws-pidfile
nodaemon=true
minfds=1024
minprocs=200
identifier=supervisor-jaws
directory=/tmp/

[supervisorctl]
prompt=supervisor-jaws
serverurl=http://localhost:9001

[program:jaws-site-daemon]
command=jaws-site --log /tmp/site-daemon.log --config /etc/jaws-site.conf --log-level DEBUG daemon
numprocs=1
startsecs=3
startretries=3
autorestart=unexpected

[program:jaws-site-jtm-rpc]
command=jaws-site --log /tmp/site-jtm-rpc.log --config /etc/jaws-site.conf --log-level DEBUG jtm-rpc
numprocs=1
startsecs=3
startretries=3
autorestart=unexpected

[program:jaws-site-central-rpc]
command=jaws-site --log /tmp/site-central-rpc.log --config /etc/jaws-site.conf --log-level DEBUG central-rpc
numprocs=1
startsecs=3
startretries=3
autorestart=unexpected

$SUPERVISOR_JAWS_CONF
EOM

exec $1

