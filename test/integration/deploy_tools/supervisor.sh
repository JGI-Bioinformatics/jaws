
function stop_services {
    $JAWS_BIN_DIR/supervisorctl stop "jaws-site:*" || true
}

function start_services {
    $JAWS_BIN_DIR/supervisorctl start "jaws-site:*" || true
}


function supervisor {
  local strategy="${1:-}"
  case "$strategy" in
    "start")
      start_services
      ;;
    "stop")
      stop_services
      ;;
    *)
      echo "Unknown strategy: $strategy"
      exit 1
      ;;
  esac
}

supervisor "$@"