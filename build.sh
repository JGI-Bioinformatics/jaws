#!/bin/bash

build_image() {
    if [ -z $1 ] && [ -z $2 ];
    then
        echo "Usage: $0 <jaws_service> <version> <env> (<uid> <gid>)"
        exit
    fi

    service_name=$1
    version=$2
    env=$3

    if [ ! -z $4 ] ;
    then
        JAWS_UID=$4
    fi

    if [ ! -z $4 ];
    then
        JAWS_GID=$5
    fi

    if [ ! -z $JAWS_UID ];
    then
      BUILD_ARGS="${BUILD_ARGS} --build-arg JAWS_UID=$JAWS_UID"
    fi

    if [ ! -z $JAWS_GID ];
    then
      BUILD_ARGS="${BUILD_ARGS} --build-arg JAWS_GID=${JAWS_GID:-75388}"
    fi

    BUILD_ARGS="--build-arg JAWS_VERSION=${version}-${env}"
    TARGET="--target ${env}"
    DIR=$service_name

    cd $DIR
    echo "****************************"
    echo "BUILDING image_version"
    echo "****************************"
    git log -n 1 --pretty="commit_count:  $(git rev-list HEAD --count)%ncommit_hash:   %h%nsubject:       %s%ncommitter:     %cN <%ce>%ncommiter_date: %ci%nauthor:        %aN <%ae>%nauthor_date:   %ai%nref_names:     %D" > image_version.yml
    cat image_version.yml

    docker build -t $service_name:$version-$env $BUILD_ARGS $TARGET .
}

build_image $1 $2 $3 $4

