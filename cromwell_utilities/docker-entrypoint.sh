#!/bin/bash

cat << EOM > /etc/cromwell.conf
include required(classpath("application"))
webservice
{
  port = $CROMWELL_PORT
}
system
{
  abort-jobs-on-terminate = true
  graceful-server-shutdown = true
  workflow-restart = false
  max-concurrent-workflows = 100000
  max-workflow-launch-count = 100000
  new-workflow-poll-rate = 1
  number-of-workflow-log-copy-workers = 20
  number-of-cache-read-workers = 50
  job-rate-control
  {
    jobs = 1
    per = 1 second
  }
}
workflow-options
{
  workflow-log-dir: "$CROMWELL_WORKFLOW_LOGS_DIR"
  workflow-log-temporary: false
  workflow-failure-mode: "ContinueWhilePossible"
  default
  {
    workflow-type: WDL
    workflow-type-version: "draft-2"
  }
}
call-caching
{
  enabled = false
  invalidate-bad-cache-result = true
}
# this is required for shifter to find image from its registry.
docker {
    hash-lookup {
        enabled = false
    }
}
backend
{
  default = "JTM"
  providers
  {
    JTM
    {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        exit-code-timeout-seconds = 60
        runtime-attributes = """
        String? docker
        String time = "00:30:00"
        Int cpu = 1
        String mem = "5G"
        String cluster = "${JAWS_SITE,,}"
        String poolname = "small"
        String constraint = "$SITE_CLUSTER_CONSTRAINT"
        String qos = "$SITE_CLUSTER_QOS"
        String account = "$SITE_CLUSTER_ACCOUNT"
        Int node = 1
        Int nwpn = 1
        Int shared = 0
        """

        kill = "jtm --config=/etc/jaws-jtm.conf kill -tid \${job_id}"

        check-alive = "jtm --config=/etc/jaws-jtm.conf isalive -tid \${job_id}"

        job-id-regex = "JTM task ID (\\\d+)"

        submit = """
CONSTRAINT=\${constraint}
[ -z \$CONSTRAINT ] || CONSTRAINT="-C \$CONSTRAINT"
jtm --config=/etc/jaws-jtm.conf \
submit \
-cmd '/bin/bash \${script}' \
-cl \${cluster} \
-t \${time} \
-c \${cpu} \
-m \${mem} \
-p \${poolname} \
-N \${node} \
-nwpn \${nwpn} \
-jid \${job_name} \
--qos \${qos} \
-A \${account} \
--shared \${shared} \
\$CONSTRAINT
"""


        submit-docker = """
CONSTRAINT=\${constraint}
$CONTAINER_CHECK
[ -z \$CONSTRAINT ] || CONSTRAINT="-C \$CONSTRAINT"
jtm --config=/etc/jaws-jtm.conf \
submit \
-cmd '$CONTAINER_WRAPPER \${docker} \${job_shell} \${script}' \
-cl \${cluster} \
-t \${time} \
-c \${cpu} \
-m \${mem} \
-p \${poolname} \
-N \${node} \
-nwpn \${nwpn} \
-jid \${job_name} \
--qos \${qos} \
-A \${account} \
--shared \${shared} \
\$CONSTRAINT
"""

        # Root directory where Cromwell writes job results in the container. This value
        # can be used to specify where the execution folder is mounted in the container.
        # it is used for the construction of the docker_cwd string in the submit-docker
        # value above AND in the generation of the "script" file.
        dockerRoot = $CROMWELL_EXECUTIONS_DIR
      }
    }
  }
}
database
{
  profile = "slick.jdbc.MySQLProfile$"
  db
  {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://$SITE_DB_HOST:$SITE_DB_PORT/cromwell_${JAWS_SITE,,}_$DEPLOYMENT_NAME?rewriteBatchedStatements=true&useSSL=false&autoReconnect=true&useUnicode=true&useJDBCCompliantTimezoneShift=true&useLegacyDatetimeCode=false&serverTimezone=UTC&allowPublicKeyRetrieval=true"
    user = "jaws"
    password = "$SITE_DB_PW"
    connectionTimeout = 5000
  }
  insert-batch-size = 2000
}
EOM

cat -n /etc/cromwell.conf
/bin/bash -c "java $JAVA_OPTS -jar /app/cromwell.jar $1"
