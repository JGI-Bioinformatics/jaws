FROM python:3.8-buster as base
ARG JAWS_UID=75388
ARG JAWS_GID=75388
RUN groupadd -g ${JAWS_GID} jaws && useradd -u ${JAWS_UID} -g ${JAWS_GID} -c  "JAWS User" --no-create-home jaws

FROM base as builder
ARG JAWS_VERSION
ENV JAWS_VERSION=${JAWS_VERSION}
RUN pip install wheel

WORKDIR /usr/rpc
COPY requirements.txt .
RUN  pip install -r requirements.txt
COPY . .

FROM builder as dev
RUN pip install -e .
WORKDIR /

FROM builder as prod
RUN python setup.py bdist_wheel
RUN pip install dist/*
WORKDIR /
