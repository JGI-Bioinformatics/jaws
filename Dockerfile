FROM python:3.10 as builder
ARG JAWS_UID=75388
ARG JAWS_GID=75388
RUN apt-get update && apt-get -y install rsync build-essential
RUN groupadd -g ${JAWS_GID} jaws && useradd -u ${JAWS_UID} -g ${JAWS_GID} -c  "JAWS User" --no-create-home jaws

WORKDIR /usr/app
COPY . /usr/app/
RUN pip install --upgrade pip \
    && pip install --no-cache -r requirements.txt && \

FROM builder as test-rpc
WORKDIR /usr/app
RUN pip install .[dev]
CMD make test-rpc

FROM builder as test-site
WORKDIR /usr/app
RUN pip install .[dev]
CMD make test-site

FROM builder as site
WORKDIR /usr/app
COPY image_version.yml image_version.yml
RUN pip install -e .

ENTRYPOINT ["jaws-site", "--config", "/etc/config/site/jaws-site.conf"]
CMD ["--log", "/var/log/rpc-server.log", "--log-level", "DEBUG", "rpc-server"]
