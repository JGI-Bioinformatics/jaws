#!/bin/bash -l

#function fix_perms() {
#    local GROUP="\$1"
#    local DIR="\$2"
#    chmod -R a+rX "\$DIR"
#    chmod -R ug+rwX "\$DIR"
#    chgrp -R "\$GROUP" "\$DIR"
#    find "\$DIR" -type d -exec chmod g+s '{}' \;
#}
#
#function set_jaws_acl() {
#    local DIR="\$1"
#    chmod 2775 \$DIR
#    setfacl -m g:$JAWS_GROUP:rwx \$DIR
#    setfacl -m g:$JAWS_USERS_GROUP:rwx \$DIR
#}

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export PYTHONIOENCODING=utf-8

#test -d "$JAWS_SCRATCH_DIR" || mkdir -p "$JAWS_SCRATCH_DIR"
#chmod 2775 "$JAWS_SCRATCH_DIR"
#test -d $JAWS_INPUTS_DIR || mkdir -p $JAWS_INPUTS_DIR
#chmod 0775 $JAWS_INPUTS_DIR
#set_jaws_acl $JAWS_INPUTS_DIR

source "$JAWS_VENV_DIR/bin/activate"
exec jaws-site --log "$JAWS_LOGS_DIR/site-transfer-daemon.log" --config "$JAWS_CONFIG_DIR/jaws-site.conf" --log-level $JAWS_LOG_LEVEL transfer-daemon
