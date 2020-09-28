#!/bin/bash

build_image() {
    if [ -z $1 ] && [ -z $2 ];
    then
        echo "Usage: $0 <jaws_service> <version>"
        exit
    fi

    service_name=$1
    version=$2

    if [ $service_name != 'cromwell' ] && [ $service_name != 'rabbitmq' ]; then
        cp -R rpc/ $service_name/rpc/
    fi
    if [ $service_name = 'cromwell' ]; then
        cd 'cromwell_utilities'
    else
        cd $service_name
    fi
    docker image build -t jaws_$service_name:$version .
    rm -rf rpc/
}

build_image $1 $2
