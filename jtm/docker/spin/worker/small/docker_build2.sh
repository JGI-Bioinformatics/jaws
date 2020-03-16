#!/bin/bash
set -e
#docker build -t sulsj/jtm:$1 . && docker push sulsj/jtm:$1
#docker build --no-cache -t sulsj/jtm:$1 . && docker push sulsj/jtm:$1
#docker build -t sulsj/jtm:$1 . && docker push sulsj/jtm:$1

if [ -z $1 ] || [ -z $2 ]; then
  echo "docker_build.sh <workertype> <cacheflag>";
  exit 1;
fi

if [ $2 == 'no-cache' ]; then
  docker build --no-cache -t registry.spin.nersc.gov/sulsj/jtm-worker-$1 - < Dockerfile.$1\
  && docker image push registry.spin.nersc.gov/sulsj/jtm-worker-$1;
else
  docker build -t registry.spin.nersc.gov/sulsj/jtm-worker-$1 - < Dockerfile.$1 \
  && docker image push registry.spin.nersc.gov/sulsj/jtm-worker-$1;
fi
