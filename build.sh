#!/bin/sh

# Build script used to help build jaws-site. It'll take care of build-args and also
# generates a image_version.yml file. Otherwise it will use the default one located
# in the repository.

local=0


help() {
    echo -n "Usage ${0} [options] [docker|apptainer] <tag>

    This is a build script for building the docker image for
    JAWS site. Use it for creating images locally.

    You can specify the JAWS user in the container.

    Options:
     -h, --help             Display this help and exit
     -u, --user <UID>       The user id of the images
     -g, --group <GID>      The group id of the images
    "
}

build_image() {
    for arg in "$@"; do
        shift
        case "$arg" in
            "--help")   set -- "$@" "-h" ;;
            "--user")  set -- "$@" "-u" ;;
            "--group")  set -- "$@" "-g" ;;
            *)          set -- "$@" "$arg"
        esac
    done

    while getopts "hlu:g:" opt; do
      case "$opt" in
        h) help >&2; exit ;;
        u) JAWS_UID=${OPTARG};;
        g) JAWS_GID=${OPTARG};;
        *) die "invalid option passed, run with -h for help." ;;
      esac
    done

    shift "$((OPTIND-1))"

    if [ -z $1 ] && [ -z $2 ];
    then
        help >&2; exit
    fi

    oci_runtime=$1
    tag=$2

    if [ ! -z $JAWS_UID ];
    then
      BUILD_ARGS="${BUILD_ARGS} --build-arg JAWS_UID=${JAWS_UID}"
    fi

    if [ ! -z $JAWS_GID ];
    then
      BUILD_ARGS="${BUILD_ARGS} --build-arg JAWS_GID=${JAWS_GID:-75388}"
    fi

    # CREATE image_version.yml
    echo "****************************"
    echo "BUILDING image_version"
    echo "****************************"
    git log -n 1 --pretty="commit_count:  $(git rev-list HEAD --count)%ncommit_hash:   %h%nsubject:       %s%ncommitter:     %cN <%ce>%ncommiter_date: %ci%nauthor:        %aN <%ae>%nauthor_date:   %ai%nref_names:     %D" > image_version.yml
    cat image_version.yml

    IMAGE_NAME="jaws-site:$tag"

    if [ "$oci_runtime" == "apptainer" ]; then
      apptainer build "jaws-site-$tag.sif" "jaws-site.def"
    else
      docker build -t $IMAGE_NAME $BUILD_ARGS .
    fi
}

build_image "$@"

