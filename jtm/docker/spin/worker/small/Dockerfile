FROM ubuntu:16.04
MAINTAINER Seung-Jin Sul ssul@lbl.gov

# Python & tools
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install -y python-dev python-pip build-essential python-setuptools
RUN apt-get install -y git zip wget sudo gawk curl file virtualenv bash
RUN pip install --upgrade pip

# Setup a virtualenv
COPY requirements.txt /
RUN apt-get update && apt-get install -y libmysqlclient-dev
RUN pip install MySQL-python cx_Oracle && \
    pip install -r requirements.txt && \
    rm -rf /var/cache/apk/*

# Clone jtm
#ARG SHA=LATEST
# v1.8.2
#ARG SHA=d5b43e848bdedbd26259b7780f9627807f51ad88
#RUN SHA=${SHA} \
#    git clone -b master http://gitlab+deploy-token-23168:EpTeXprf8z3zByqHyvnx@gitlab.com/sulsj/jtm_task_manager.git
RUN git clone -b master http://gitlab+deploy-token-23168:EpTeXprf8z3zByqHyvnx@gitlab.com/sulsj/jtm_task_manager.git
# RUN git clone -b test_docker http://gitlab+deploy-token-23168:EpTeXprf8z3zByqHyvnx@gitlab.com/sulsj/jtm_task_manager.git

RUN cd /usr/local/bin/ && \
    ln -s /jtm_task_manager/jgi-task-manager && \
    ln -s /jtm_task_manager/jtm-worker && \
    ln -s /jtm_task_manager/jtm-submit && \
    ln -s /jtm_task_manager/jtm-status && \
    ln -s /jtm_task_manager/jtm-isalive && \
    ln -s /jtm_task_manager/jtm-check-manager && \
    ln -s /jtm_task_manager/jtm-check-worker && \
    ln -s /jtm_task_manager/jtm-remove-pool && \
    ln -s /jtm_task_manager/jtm-kill

RUN mkdir -p /logs
RUN mkdir -p /jobs

# EXPOSE 4369 5671 5672 15671 15672 25672 55672 3306

# COPY run_jtm-worker.sh /usr/local/bin/
# COPY run_jtm.sh /usr/local/bin/
ENTRYPOINT ["python", "/jtm_task_manager/jtm-worker", "-tp", "small", "-jd", "/jobs", "-ld", "/logs"] # configures a container that will run as an executable.
# ENTRYPOINT ["python", "/jtm_task_manager/jtm-worker", "-tp", "short"] # configures a container that will run as an executable.
# CMD ["jtm-worker"] # sets default command and/or parameters
