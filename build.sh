#!/bin/sh

local=0


help() {
    echo -n "Usage ${0} [options] <service-name> <version> <environment>

    This is a build script for building the docker images for
    JAWS services. Use it for creating images locally or in a
    CI pipeline to automatically create and upload container
    images to Gitlab Repository.

    Takes in a service name (central, site), a version, and the environment
    (dev, prod). A development environment will create a development container.

    Options:
     -h, --help             Display this help and exit
     -l, --local            Build local images without a registry
     -u, --user <UID>       The user id of the images
     -g, --group <GID>      The group id of the images
    "
}

build_image() {
    for arg in "$@"; do
        shift
        case "$arg" in
            "--help")   set -- "$@" "-h" ;;
            "--local")  set -- "$@" "-l" ;;
            "--user")  set -- "$@" "-u" ;;
            "--group")  set -- "$@" "-g" ;;
            *)          set -- "$@" "$arg"
        esac
    done

    while getopts "hlu:g:" opt; do
      case "$opt" in
        h) help >&2; exit ;;
        l) local=1;;
        u) JAWS_UID=${OPTARG};;
        g) JAWS_GID=${OPTARG};;
        *) die "invalid option passed, run with -h for help." ;;
      esac
    done

    shift "$((OPTIND-1))"

    if [ -z $1 ] && [ -z $2 ] && [ -z $3 ];
    then
        help >&2; exit
    fi

    service_name=$1
    version=$2
    env=$3

    GITLAB_REGISTRY="code.jgi.doe.gov:5050/advanced-analysis/jaws-site"

    if [ ! -z $JAWS_UID ];
    then
      BUILD_ARGS="${BUILD_ARGS} --build-arg JAWS_UID=${JAWS_UID}"
    fi

    if [ ! -z $JAWS_GID ];
    then
      BUILD_ARGS="${BUILD_ARGS} --build-arg JAWS_GID=${JAWS_GID:-75388}"
    fi

    BUILD_ARGS="--build-arg JAWS_VERSION=${version}"
    TARGET="--target ${env}"
    DIR=$service_name

    cd $DIR

    IMAGE_NAME="$service_name:$version"
    if [ $local -eq 0 ];
    then
        IMAGE_NAME=$GITLAB_REGISTRY/$service_name:$version
    fi

    docker build -t $IMAGE_NAME $BUILD_ARGS $TARGET .

    if [ $local -eq 0 ];
    then
        docker image push $IMAGE_NAME
    fi
}

build_image $@

