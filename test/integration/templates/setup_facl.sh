#!/bin/bash

# Creates the inputs directory that sets the ACL for jaws service and user groups.
# Make sure to run this script prior to deployment so that ACL rules are properly set.

echo "Setting up jaws-$JAWS_DEPLOYMENT_NAME"

DIR="$JAWS_SCRATCH_DIR/jaws-$JAWS_DEPLOYMENT_NAME"
test -d "$DIR" || mkdir "$DIR"
chgrp "$JAWS_GROUP" "$DIR"
chmod 0775 "$DIR"
setfacl -m d:g:"$JAWS_GROUP":rwx "$DIR"
setfacl -m d:g:"$JAWS_USERS_GROUP":rwx "$DIR"

DIR="$JAWS_SCRATCH_DIR/jaws-$JAWS_DEPLOYMENT_NAME/cromwell-executions"
test -d "$DIR" || mkdir "$DIR"
chgrp "$JAWS_GROUP" "$DIR"
chmod 0775 "$DIR"
setfacl -m d:g:"$JAWS_GROUP":rwx "$DIR"
setfacl -m d:g:"$JAWS_USERS_GROUP":rx "$DIR"

DIR="$JAWS_SCRATCH_DIR/jaws-$JAWS_DEPLOYMENT_NAME/inputs"
test -d "$DIR" || mkdir "$DIR"
chgrp "$JAWS_GROUP" "$DIR"
chmod 0775 "$DIR"
setfacl -m d:g:"$JAWS_GROUP":rwx "$DIR"
setfacl -m d:g:"$JAWS_USERS_GROUP":rwx "$DIR"

DIR="$JAWS_SCRATCH_DIR/jaws-$JAWS_DEPLOYMENT_NAME/downloads"
test -d "$DIR" || mkdir "$DIR"
chgrp "$JAWS_GROUP" "$DIR"
chmod 0775 "$DIR"
setfacl -m d:g:"$JAWS_GROUP":rwx "$DIR"
setfacl -m d:g:"$JAWS_USERS_GROUP":rx "$DIR"
