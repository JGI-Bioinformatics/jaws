#!/bin/bash
set -e
#docker build -t sulsj/jtm:$1 . && docker push sulsj/jtm:$1
#docker build --no-cache -t sulsj/jtm:$1 . && docker push sulsj/jtm:$1
#docker build -t sulsj/jtm:$1 . && docker push sulsj/jtm:$1

docker build --no-cache -t registry.spin.nersc.gov/sulsj/jtm-worker-small . \
&& docker image push registry.spin.nersc.gov/sulsj/jtm-worker-small
#docker build -t registry.spin.nersc.gov/sulsj/jtm-worker-small . \
#&& docker image push registry.spin.nersc.gov/sulsj/jtm-worker-small
