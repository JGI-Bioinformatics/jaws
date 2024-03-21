#!/bin/bash -l

${JAWS_APPTAINER_PATH} run --cleanenv --env-file="$JAWS_CONFIG_DIR/site.env" --no-home --mount src="${JAWS_CONFIG_DIR}/jaws-site.conf",dst=/etc/config/site/jaws-site.conf,ro --mount src="${JAWS_LOGS_DIR}",dst=/var/log/site "${JAWS_BIN_DIR}/site-${JAWS_SITE_VERSION}.sif" --log "/var/log/site/${SERVICE}.log" --log-level "${JAWS_LOG_LEVEL}" "${SERVICE}"

