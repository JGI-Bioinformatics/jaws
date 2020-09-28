FROM giovtorres/docker-centos7-slurm:latest


COPY rpc/ rpc/

COPY docker-entrypoint.sh bin/


RUN pip3.7 install wheel && cd rpc/ && python3.7 setup.py bdist_wheel && pip3.7 install dist/*

WORKDIR jtm/
COPY jaws_jtm/ jaws_jtm/
COPY tests/ tests/
COPY requirements.txt .
COPY setup.py .
RUN python3.7 -m venv --system-site-packages jaws_jtm && python3.7 -m venv --system-site-packages jaws_jtm_worker
RUN source jaws_jtm/bin/activate && pip install wheel && python3.7 setup.py bdist_wheel && pip3.7 install dist/* && deactivate
RUN source jaws_jtm_worker/bin/activate && python3.7 setup.py bdist_wheel && pip3.7 install dist/* && deactivate

ENTRYPOINT ["docker-entrypoint.sh"]

CMD ["jtm --config=/jtm/jaws_jtm/jaws-jtm.conf --debug manager -ld /tmp/"]
