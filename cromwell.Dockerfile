FROM ubuntu:focal

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -q && \
    apt-get upgrade -qq && \
    apt-get install -y curl wget gnupg gnupg2 git openjdk-11-jdk && \
    rm -rf /var/lib/apt/lists/*

# docker build -t vanessa/cromwell-dev .

# Development environment for Cromwell that includes:
#
#   Scala 2.13
#   SBT 1.x
#   Java 11
#   Git

# Env variables
ENV SCALA_VERSION 2.13.8
ENV SBT_VERSION 1.5.5
ENV SBT_OPTS="-Xmx2G"

#
## Scala
#

RUN mkdir -p /home/pigman && \
    curl \
    --location --fail --silent --show-error \
    https://downloads.typesafe.com/scala/$SCALA_VERSION/scala-$SCALA_VERSION.tgz | tar xfz - -C /opt/ && \
    echo >> /home/pigman/.bashrc && \
    echo "export PATH=/opt/scala-$SCALA_VERSION/bin:$PATH" >> /home/pigman/.bashrc

#
## sbt
#

# non-deb package installation instructions adapted from
# - https://github.com/sbt/sbt/releases/tag/v1.4.9
# - https://github.com/broadinstitute/scala-baseimage/pull/4/files
RUN curl \
    --location --fail --silent --show-error \
    "https://github.com/sbt/sbt/releases/download/v$SBT_VERSION/sbt-$SBT_VERSION.tgz" | \
    tar zxf - -C /usr/share && \
    update-alternatives --install /usr/bin/sbt sbt /usr/share/sbt/bin/sbt 1 && \
    sbt -Dsbt.supershell=false -Dsbt.rootdir=true sbtVersion

# Instruct user to add code here during development
RUN useradd pigman && \
    mkdir -p /code
WORKDIR /code

RUN git clone --branch 86 --depth 1 https://github.com/broadinstitute/cromwell.git
WORKDIR /code/cromwell
RUN sbt assembly

FROM debian:bookworm

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -q && \
    apt-get upgrade -qq && \
    apt-get -qq -y install podman && \
    apt-get install -y curl wget fuse3 fuse-overlayfs libcap2-bin vim openjdk-17-jdk containers-storage iptables && \
    rm -rf /var/lib/apt/lists/*

RUN useradd -m -s /bin/bash -d /home/cromwell cromwell; \
    echo cromwell:10000:5000 > /etc/subuid; \
    echo cromwell:10000:5000 > /etc/subgid;

VOLUME /var/lib/containers
VOLUME /home/podman/.local/share/containers

ADD https://raw.githubusercontent.com/containers/libpod/master/contrib/podmanimage/stable/containers.conf /etc/containers/containers.conf
ADD https://raw.githubusercontent.com/containers/libpod/master/contrib/podmanimage/stable/podman-containers.conf /home/podman/.config/containers/containers.conf

RUN ls /etc/containers/
RUN ln -s /usr/bin/podman /usr/bin/docker

RUN mkdir code
RUN chown cromwell:cromwell -R /home/cromwell
RUN chown -R cromwell:cromwell -R /code

RUN cp /usr/share/containers/storage.conf /etc/containers/
RUN chmod 644 /etc/containers/containers.conf; sed -i -e 's|^#mount_program|mount_program|g' -e '/additionalimage.*/a "/var/lib/shared",' -e 's|^mountopt[[:space:]]*=.*$|mountopt = "nodev,fsync=0"|g' /etc/containers/storage.conf
RUN mkdir -p /var/lib/shared/overlay-images /var/lib/shared/overlay-layers /var/lib/shared/vfs-images /var/lib/shared/vfs-layers; touch /var/lib/shared/overlay-images/images.lock; touch /var/lib/shared/overlay-layers/layers.lock; touch /var/lib/shared/vfs-images/images.lock; touch /var/lib/shared/vfs-layers/layers.lock

ENV _CONTAINERS_USERNS_CONFIGURED=""

COPY --from=0 /code/cromwell/server/target/scala-2.13/cromwell-86*.jar /code/cromwell.jar
WORKDIR code
CMD ["java", "-jar", "/code/cromwell.jar", "server"]
