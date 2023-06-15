#!/bin/bash

# This script is only used when jaws-site is submitted to a cluster scheduler (e.g. SLURM)
# rather than being installed on a dedicated server.

JAWS_BIN_DIR="$JAWS_BIN_DIR"
JAWS_LOGS_DIR="$JAWS_LOGS_DIR"
JAWS_GITLAB_RUNNER="$JAWS_GITLAB_RUNNER"
JAWS_GITLAB_RUNNER_CONFIG="$JAWS_GITLAB_RUNNER_CONFIG"
LOGFILE="$JAWS_LOGS_DIR/jaws-site.log"


# Dynamic dns setting
HOSTIP=`hostname -i`
SITEDNSNAME="${JAWS_SITE_DNS_NAME}"
PSCRATCHDNSNAME="pscratch-s3.jaws.dyn.jgi.doe.gov"

# The last Minio release that allows the creation of new NAS gateways was docker pull 2022-05-26T05-48-41Z
# The last release that allowed you to use the NAS gateway on an existing setup was 2022-10-24T18-35-07Z
# Export these so that Minio sees them
#####################################   Steve Chan wants this here uncommented:
#export MINIO_ROOT_USER=jaws
#export MINIO_ROOT_PASSWORD=thisisthejawspassword.
#export MINIO_BROWSER=off

#MINIO_EXE="/global/u2/j/jaws/minio/minio"
#MINIOPORT=9000
#LISTEN="${DOLLAR}{PSCRATCHDNSNAME}:${DOLLAR}{MINIOPORT}"

# See if we have a HTTP listener on the existing DNSNAME port 9000, if so then bail since that
# needs to be shutdown first
#curl http://${DOLLAR}{LISTEN}/
#if [ ${DOLLAR}? -eq 0 ]; then
    #echo "Existing listener at ${DOLLAR}{LISTEN}. It must be shutdown before we start another"
    #pkill -u jaws minio
#fi
######################################

echo "Setting DNS for ${DOLLAR}SITEDNSNAME with ${DOLLAR}HOSTIP"
TSIGPATH=${TSIGPATH}
cat <<EOF >/tmp/JAWS_SITEDNS.${DOLLAR}${DOLLAR}
server nsx.lbl.gov
zone dyn.jgi.doe.gov
update del ${DOLLAR}SITEDNSNAME A
update add ${DOLLAR}SITEDNSNAME 30 A ${DOLLAR}HOSTIP
send
EOF

nsupdate -v -k ${DOLLAR}{TSIGPATH} </tmp/JAWS_SITEDNS.${DOLLAR}${DOLLAR}
if [ ${DOLLAR}? -eq 0 ]; then
    sleep 30
    echo "Successfully set DNS entry"
    dig +short ${DOLLAR}SITEDNSNAME
    #${DOLLAR}MINIO_EXE gateway nas --address ${DOLLAR}HOSTIP:${DOLLAR}MINIOPORT /pscratch/sd/j/jaws
else
    echo "Failed to set DNS entry. Exiting"
    exit 1
fi


echo "Starting jaws-site on `hostname -i`" > "${DOLLAR}LOGFILE"

# Start supervisord
"${DOLLAR}JAWS_BIN_DIR/supervisord" || true
"${DOLLAR}JAWS_BIN_DIR/supervisorctl" start all || true
