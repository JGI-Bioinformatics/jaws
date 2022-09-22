FROM python:3.8-buster as builder
ARG JAWS_UID=75388
ARG JAWS_GID=75388
RUN groupadd -g ${JAWS_GID} jaws && useradd -u ${JAWS_UID} -g ${JAWS_GID} -c  "JAWS User" --no-create-home jaws

WORKDIR /usr/app
COPY rpc rpc
RUN cd rpc && python setup.py install
COPY site/requirements.txt requirements.txt
RUN pip install -r requirements.txt

FROM builder as site
WORKDIR /usr/app
COPY site .
RUN python setup.py install
USER jaws