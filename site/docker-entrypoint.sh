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
password = $RMQ_PW
host = $RMQ_HOST
port = $RMQ_PORT
vhost = jaws_$DEPLOYMENT_NAME
queue = JAWS_$JAWS_SITE
num_threads = 5
max_retries = 3
[CENTRAL_RPC_CLIENT]
user = jaws
password = $RMQ_PW
host = $RMQ_HOST
port = $RMQ_PORT
vhost = jaws_$DEPLOYMENT_NAME
queue = CENTRAL
[GLOBUS]
client_id = $GLOBUS_CLIENT_ID
client_secret = $GLOBUS_CLIENT_SECRET
host_path = $GLOBUS_HOST_PATH
endpoint_id = $GLOBUS_EP
[SITE]
id = $JAWS_SITE
uploads_dir = $UPLOADS_DIR
[CROMWELL]
url = http://cromwell:$CROMWELL_PORT
EOM
exec $1

