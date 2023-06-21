#!/bin/bash -l

function login_gitlab_registry {
  local registry="$1"
  echo "$CI_REGISTRY_PASSWORD" | apptainer remote login --username "$CI_REGISTRY_USER" --password-stdin "$registry"
}

function pull_container {
  local registry="$1"
  local image="$2"
  local tag="$3"
  echo "Pulling container"
  if [ $CONTAINER_RUNTIME == "apptainer" ]; then
    apptainer-pull --force "${JAWS_BIN_DIR}/site-${JAWS_SITE_VERSION}.sif" "$registry/$image:$tag"
  else
    echo "Unknown container runtime: $CONTAINER_RUNTIME"
    exit 1
  fi
}

function install_venv {
  local python="$1"
  local venv_dir="$2"
  test -d  $venv_dir && rm -rf $venv_dir
  [[ -n "$JAWS_LOAD_PYTHON" ]] && $JAWS_LOAD_PYTHON
  $python -m venv $venv_dir && \
  . "$venv_dir/bin/activate" && \
  pip install build && \
  pip install wheel && \
  make pkg && \
  pip install rpc/dist/*.whl && \
  pip install site/dist/*.whl && \
  deactivate
}

function install_jaws_site {
  local install_method="$1"
  case "$install_method" in
    "venv")
      install_venv "$JAWS_PYTHON" "$JAWS_VENV_DIR"
      ;;
    "apptainer")
      login_gitlab_registry "oras://$CI_REGISTRY"
      pull_container "oras://$CI_REGISTRY" "advanced-analysis/jaws-site" "$JAWS_SITE_VERSION"
      ;;
    *)
      echo "Unknown install method: $install_method"
      exit 1
      ;;
  esac
}