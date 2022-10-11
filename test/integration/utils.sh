# Helper functions to be used in deployment

function setup_dir() {
    local DIR="$1"
    local GROUP="$2"
    local PERMS="$3"
    [[ -z $DIR ]] && echo "missing DIR" && exit 1
    [[ -z $GROUP ]] && echo "missing GROUP" && exit 1
    echo "... mkdir $DIR"
    test -d "$DIR" || mkdir -p "$DIR"
    chgrp "$GROUP" "$DIR"
    chmod "$PERMS" "$DIR"
}


function setup_dirs() {
    local FOLDERS=$1
    local GROUP=$2
    local PERMS=$3
    for DIR in $FOLDERS; do
      setup_dir "${!DIR}" "$GROUP" "$PERMS"
    done
}


# If any of the specified variables are undefined,
# print all undefined variable names to stderr and exit,
# otherwise return silently.
function validate_vars {
    local REQUIRED_VARS=$1
    local RESULT=0
    for VAR in $REQUIRED_VARS; do
        if [[ -z "${!VAR+xxx}" ]]; then
            echo "Env var not defined: $VAR">&2
            RESULT=1
        fi
    done
    if [ $RESULT -eq 1 ]; then
        exit 1
    fi
}


function fix_perms() {
    local GROUP=$1
    local DIR="$2"
    local GROUP_PERMS=${3:-rwX}
    chmod -R a+rX "$DIR"
    chmod -R u+rwX "$DIR"
    chmod -R g=$GROUP_PERMS "$DIR"
    chgrp -R $GROUP "$DIR"
    find "$DIR" -type d -exec chmod g+s '{}' \;
}


function not_available() {
  command -v $1 >/dev/null 2>&1
}
