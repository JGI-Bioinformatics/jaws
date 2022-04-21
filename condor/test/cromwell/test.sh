PATH="/global/cfs/projectdirs/jaws/condor/condor/bin":"/global/cfs/projectdirs/jaws/condor/condor/sbin":$PATH
out=$(shifterimg lookup jfroula/img-omics@sha256:d90f7718861adcb81db8664e4b04452381628bdb78d0613aa959ce924d32b2c3 || shifterimg pull jfroula/img-omics@sha256:d90f7718861adcb81db8664e4b04452381628bdb78d0613aa959ce924d32b2c3)
ret=$?
if [[ $ret != 0 || $(echo $out | grep "FAILURE") ]]; then
  echo "Invalid container name or failed to pull container, jfroula/img-omics@sha256:d90f7718861adcb81db8664e4b04452381628bdb78d0613aa959ce924d32b2c3!" >&2
  exit $ret
else
  echo "Successfully pulled jfroula/img-omics@sha256:d90f7718861adcb81db8664e4b04452381628bdb78d0613aa959ce924d32b2c3!"
fi
chmod 755 /global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/script
#FLAVOR=$(/global/cfs/cdirs/m3408/aim2/cromwell2/get_flavor.sh 1 512.0 )
FLAVOR=dynamic
printenv > env.debug
#cat > /global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/dockerScript <<EOF
##!/bin/bash
##PATH=/global/common/software/m3408/cromwell:$PATH
#PATH="/global/cfs/projectdirs/jaws/condor/condor/bin":"/global/cfs/projectdirs/jaws/condor/condor/sbin":$PATH
##cd /global/cscratch1/sd/nmdcda/
##/global/cfs/cdirs/m3408/aim2/cromwell/shifter_exec.sh jfroula/img-omics@sha256:d90f7718861adcb81db8664e4b04452381628bdb78d0613aa959ce924d32b2c3 /bin/bash /global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/script
#/global/cfs/projectdirs/jaws/jaws-install/jaws-prod/shifter_exec.sh jfroula/img-omics@sha256:d90f7718861adcb81db8664e4b04452381628bdb78d0613aa959ce924d32b2c3 /global/dna/shared/databases/jaws/refdata /refdata /bin/bash /global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/script
#EOF
#chmod 755 /global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/dockerScript
#cat > /global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/submitFile <<EOF
#Iwd=/global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution
#+Owner=UNDEFINED
##requirements=$FLAVOR
##leave_in_queue=true
##request_memory=512.0
##request_disk=25600.0
##request_cpus=1
##priority=
#error=/global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/stderr
#output=/global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/stdout
#log_xml=true
#executable=/global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/dockerScript
#log=/global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/execution.log
#queue
#EOF
#condor_submit /global/homes/j/jaws_jtm/ssul/condor/cromwell/cromwell-executions/fq_count/68f73673-f699-4980-801b-5fb3adccb18c/call-count_seqs/execution/submitFile
