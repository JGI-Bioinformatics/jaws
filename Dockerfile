FROM python:3.11 as builder
ARG JAWS_UID=75388
ARG JAWS_GID=75388
RUN apt-get update && apt-get -y install rsync
RUN groupadd -g ${JAWS_GID} jaws && useradd -u ${JAWS_UID} -g ${JAWS_GID} -c  "JAWS User" --no-create-home jaws

WORKDIR /usr/app
COPY rpc rpc
RUN cd rpc && pip install --upgrade pip  \
    && pip install -r requirements.txt && \
    pip install .
COPY site/requirements.txt requirements.txt
RUN pip install -r requirements.txt

FROM builder as test-rpc
WORKDIR /usr/app/rpc
RUN pip install -r dev-requirements.txt
CMD python -m pytest --cov=jaws_rpc --junitxml=coverage/rpc.xml tests/ && coverage xml

FROM builder as test-site
WORKDIR /usr/app
COPY site .
RUN pip install -r dev-requirements.txt
CMD python -m pytest --cov=jaws_site --junitxml=coverage/site.xml tests/ && coverage xml

FROM builder as site
WORKDIR /usr/app
COPY site .
RUN pip install .

ENTRYPOINT ["jaws-site", "--config", "/etc/config/site/jaws-site.conf"]
CMD ["--log", "/var/log/rpc-server.log", "--log-level", "DEBUG", "rpc-server"]