#!/usr/bin/env bash

# Installs required dependencies inside the docker image used for publishing cromwell

set -eou pipefail

apt update
apt install \
    apt-transport-https \
    curl \
    git \
    gnupg \
    wget \
    ca-certificates \
    unzip \
    zip \
    -y --no-install-recommends

apt update
apt install -y openjdk-11-jdk

# Install jq 1.6 to ensure --rawfile is supported
curl -L https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 -o /usr/bin/jq
chmod +x /usr/bin/jq

# sbt launcher non-deb package installation instructions adapted from
# - https://github.com/sbt/sbt/releases/tag/v1.4.9
# - https://github.com/broadinstitute/scala-baseimage/pull/4/files
curl --location --fail --silent --show-error "https://github.com/sbt/sbt/releases/download/v1.8.2/sbt-1.8.2.tgz" |
    tar zxf - -C /usr/share
update-alternatives --install /usr/bin/sbt sbt /usr/share/sbt/bin/sbt 1

# Update sbt launcher
sbt -Dsbt.supershell=false -Dsbt.rootdir=true sbtVersion
