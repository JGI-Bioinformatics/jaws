#!/bin/bash
apptainer instance start --cleanenv --env-file="$JAWS_CONFIG_DIR/site.env" --no-home --mount src="${JAWS_CONFIG_DIR}/jaws-site.conf",dst=/etc/config/site/jaws-site.conf,ro --mount src="${JAWS_LOGS_DIR}",dst=/var/log/site --pid-file ${JAWS_LOGS_DIR}/jaws-${SERVICE}-${JAWS_DEPLOYMENT_NAME}.pid "${JAWS_BIN_DIR}/site-${JAWS_SITE_VERSION}.sif" jaws-${SERVICE}-${JAWS_DEPLOYMENT_NAME} --log "/var/log/site/${SERVICE}.log" --log-level "${JAWS_LOG_LEVEL}" "${SERVICE}"

