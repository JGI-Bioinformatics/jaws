#!/usr/bin/env bash

CONF=""

CONF="$CONF
[SITE:LOCAL]
user = jaws
password = $JAWS_CENTRAL_RMQ_PW
host = $JAWS_CENTRAL_RMQ_HOST
port = $JAWS_CENTRAL_RMQ_PORT
vhost = jaws_$DEPLOYMENT_NAME
queue = JAWS_LOCAL
globus_endpoint = /
globus_basepath = /
uploads_subdir = /uploads
max_ram_gb = 100G
"
cat << EOM > /etc/jaws-central.conf
[JAWS]
name = $DEPLOYMENT_NAME
version = $JAWS_VERSION
docs_url = $JAWS_DOCS_URL
[DB]
dialect = mysql+mysqlconnector
host = $JAWS_CENTRAL_DB_HOST
port = $JAWS_CENTRAL_DB_PORT
user = jaws
password = $JAWS_CENTRAL_DB_PW
db = jaws_central_$DEPLOYMENT_NAME
[RPC_SERVER]
user = jaws
password = $JAWS_CENTRAL_RMQ_PW
host = $JAWS_CENTRAL_RMQ_HOST
port = $JAWS_CENTRAL_RMQ_PORT
vhost = jaws_$DEPLOYMENT_NAME
queue = CENTRAL
num_threads = 5
max_retries = 3
[GLOBUS]
client_id = $JAWS_GLOBUS_CLIENT_ID
[HTTP]
auth_port = $JAWS_AUTH_PORT
rest_port = $JAWS_REST_PORT
$CONF
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

[program:jaws-central-rest]
command=jaws-central --log /tmp/central-rest.log --config /etc/jaws-central.conf --log-level DEBUG rest
numprocs=1
startsecs=3
startretries=3
autorestart=unexpected

[program:jaws-central-auth]
command=jaws-central --log /tmp/central-auth.log --config /etc/jaws-central.conf --log-level DEBUG auth
numprocs=1
startsecs=3
startretries=3
autorestart=unexpected

[program:jaws-central-rpc]
command=jaws-central --log /tmp/central-rpc.log --config /etc/jaws-central.conf --log-level DEBUG rpc
numprocs=1
startsecs=3
startretries=3
autorestart=unexpected
$SUPERVISOR_JAWS_CONF
EOM

exec $1
