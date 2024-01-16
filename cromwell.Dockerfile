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
    apt-get -qq -y install podman
    apt-get install -y curl wget openjdk-17-jdk iptables && \
    rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/podman /usr/bin/docker
RUN mkdir code
WORKDIR code
COPY --from=0 /code/cromwell/server/target/scala-2.13/cromwell-86*.jar /code/cromwell.jar
CMD ["java", "-jar", "/code/cromwell.jar", "server"]
