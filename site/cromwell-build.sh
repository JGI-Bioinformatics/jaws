#!/usr/bin/env bash

# Builds docker image for publishing cromwell

set -euo pipefail

export DOCKER_CLI_EXPERIMENTAL=enabled

docker buildx rm cromwell-multi-arch-builder || true
docker buildx create --use --name cromwell-multi-arch-builder

git clone https://github.com/broadinstitute/cromwell.git
cp cromwell.Dockerfile cromwell/Dockerfile
cp cromwell-setup.sh cromwell/docker-setup.sh
cd cromwell

build_root="$( dirname "${BASH_SOURCE[0]}" )"
docker buildx build "${build_root}" --platform linux/arm64 -t jgi/cromwell-jaws --load 

docker buildx rm cromwell-multi-arch-builder
