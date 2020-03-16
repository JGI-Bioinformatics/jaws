#!/bin/bash
#docker build -t sulsj/jtm:$1 . && docker push sulsj/jtm:$1
docker build --no-cache -t sulsj/jtm:$1 . && docker push sulsj/jtm:$1
