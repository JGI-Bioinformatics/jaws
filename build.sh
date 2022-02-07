#!/bin/bash

set -x

build_image() {
    if [ -z $1 ] && [ -z $2 ];
    then
        echo "Usage: $0 <jaws_service> <version>"
        exit
    fi

    service_name=$1
    version=$2

    if [ $service_name = 'cromwell' ]; then
        cd 'cromwell_utilities'
    else
        cd $service_name
    fi

    git clone https://code.jgi.doe.gov/advanced-analysis/jaws
    docker build --no-cache -t  jaws_$service_name:$version --build-arg JAWS_VERSION=${version} --build-arg VERSION=2.5.0 .
    rm -rf jaws/
}

build_image $1 $2
