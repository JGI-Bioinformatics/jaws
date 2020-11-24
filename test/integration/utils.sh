# Helper functions to be used in deployment

## CREATE DIRS AND SET GROUP AND PERMISSIONS
## - there are expected to be "jaws" and "jtm" users
## - "jaws" user belongs to "jaws" and "jtm" groups
## - "jtm" user belongs to "jtm" group and no others for security purposes
## - users don't belong to either "jaws" or "jtm" groups (protected)
function setup_dir {
  local DIR="$1"
  local GROUP="$2"
  local PERMS="$3"
  [[ -z $DIR ]] && echo "missing DIR" && exit 1
  [[ -z $GROUP ]] && echo "missing GROUP" && exit 1
  echo "... $DIR"
  test -d "$DIR" || mkdir -p "$DIR"
  chgrp "$GROUP" "$DIR"
  chmod "$PERMS" "$DIR"
}


function validate_vars {
    RESULT=0
    if [[ -z ${!VAR} ]]; then
        echo "Missing env var: $VAR">&2
        RESULT=1
    else
        echo "... $VAR : OK"
    fi
    [ $RESULT -eq 0 ] || exit $RESULT
}


function set_deployment_var {
  VAR_NAME="$1"
  DEPLOYMENT_NAME="$2"
  SRC_VAR_NAME="${DEPLOYMENT_NAME}_${VAR_NAME}"
  if [ -z ${!SRC_VAR_NAME+xxx} ]; then
    echo "Missing env var: $SRC_VAR_NAME">&2
    exit 1
  fi
  VALUE="${!SRC_VAR_NAME}"
  echo "... $VAR_NAME"
  export $VAR_NAME
  printf -v "$VAR_NAME" "%s" "$VALUE"
}


function set_site_var {
  JAWS_SITE="$1"
  VAR_NAME="$2"
  SRC_VAR_NAME="${JAWS_SITE}_${VAR_NAME}"
  if [ -z ${!SRC_VAR_NAME+xxx} ]; then
    echo "Missing env var: $SRC_VAR_NAME">&2
    exit 1
  fi
  SITE_VAR_NAME="SITE_${VAR_NAME}"
  VALUE="${!SRC_VAR_NAME}"
  echo "... $SITE_VAR_NAME"
  export $SITE_VAR_NAME
  printf -v "$SITE_VAR_NAME" "%s" "$VALUE"
}


function fix_perms() {
    local GROUP=$1
    local DIR="$2"
    chmod -R a+rX "$DIR"
    chmod -R ug+rwX "$DIR"
    chgrp -R $GROUP "$DIR"
    find "$DIR" -type d -exec chmod g+s '{}' \;
}
