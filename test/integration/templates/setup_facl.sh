#!/bin/bash

test -d "$JAWS_SCRATCH_DIR" || mkdir "$JAWS_SCRATCH_DIR"
chgrp "$JAWS_GROUP" "$JAWS_SCRATCH_DIR"
chmod 0775 "$JAWS_SCRATCH_DIR"
setfacl -m d:g:"$JAWS_GROUP":rwx "$JAWS_SCRATCH_DIR"
if [[ "$JAWS_USERS_GROUP" == "" ]]; then
    setfacl -m d:o::rwx "$JAWS_SCRATCH_DIR"
else
    setfacl -m d:g:"$JAWS_USERS_GROUP":rwx "$JAWS_SCRATCH_DIR"
fi

test -d "$JAWS_CROMWELL_EXECUTIONS_DIR" || mkdir "$JAWS_CROMWELL_EXECUTIONS_DIR"
chgrp "$JAWS_GROUP" "$JAWS_CROMWELL_EXECUTIONS_DIR"
chmod 0775 "$JAWS_CROMWELL_EXECUTIONS_DIR"
setfacl -m d:g:"$JAWS_GROUP":rwx "$JAWS_CROMWELL_EXECUTIONS_DIR"
if [[ "$JAWS_USERS_GROUP" == "" ]]; then
    setfacl -m d:o::rx "$JAWS_CROMWELL_EXECUTIONS_DIR"
else
    setfacl -m d:g:"$JAWS_USERS_GROUP":rx "$JAWS_CROMWELL_EXECUTIONS_DIR"
fi

test -d "$JAWS_INPUTS_DIR" || mkdir "$JAWS_INPUTS_DIR"
chgrp "$JAWS_GROUP" "$JAWS_INPUTS_DIR"
chmod 0775 "$JAWS_INPUTS_DIR"
setfacl -m d:g:"$JAWS_GROUP":rwx "$JAWS_INPUTS_DIR"
if [[ "$JAWS_USERS_GROUP" == "" ]]; then
    setfacl -m d:o::rwx "$JAWS_INPUTS_DIR"
else
    setfacl -m d:g:"$JAWS_USERS_GROUP":rwx "$JAWS_INPUTS_DIR"
fi

test -d "$JAWS_DOWNLOADS_DIR" || mkdir "$JAWS_DOWNLOADS_DIR"
chgrp "$JAWS_GROUP" "$JAWS_DOWNLOADS_DIR"
chmod 0775 "$JAWS_DOWNLOADS_DIR"
setfacl -m d:g:"$JAWS_GROUP":rwx "$JAWS_DOWNLOADS_DIR"
if [[ "$JAWS_USERS_GROUP" == "" ]]; then
    setfacl -m d:o::rwx "$JAWS_DOWNLOADS_DIR"
else
    setfacl -m d:g:"$JAWS_USERS_GROUP":rx "$JAWS_DOWNLOADS_DIR"
fi
