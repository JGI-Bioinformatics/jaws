FROM ubuntu:rolling

WORKDIR /cromwell-publish/
COPY docker-setup.sh ./
RUN ./docker-setup.sh