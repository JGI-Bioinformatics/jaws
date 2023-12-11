FROM python:3.10-slim-bookworm as builder

ARG JAWS_UID=75388
ARG JAWS_GID=75388
RUN export DEBIAN_FRONTEND=noninteractive && \
  apt-get update && \
  apt-get -y upgrade && \
  apt-get install -y --no-install-recommends tini procps rsync npm \
  build-essential git make zip && \
  apt-get -y clean && \
  rm -rf /var/lib/apt/lists/*
RUN groupadd -g ${JAWS_GID} jaws && useradd -u ${JAWS_UID} -g ${JAWS_GID} -c  "JAWS User" --no-create-home jaws

WORKDIR /usr/app
COPY . /usr/app/
RUN make init-dev

FROM builder as test-rpc
WORKDIR /usr/app
RUN make init-dev
CMD make test-rpc

FROM builder as test-site
WORKDIR /usr/app
RUN make init-dev
CMD make test-site

FROM builder as site
WORKDIR /usr/app
COPY image_version.yml image_version.yml
RUN make init

ENTRYPOINT ["jaws-site", "--config", "/etc/config/site/jaws-site.conf"]
CMD ["--log", "/var/log/rpc-server.log", "--log-level", "DEBUG", "rpc-server"]
