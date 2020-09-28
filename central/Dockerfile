FROM python:3.8-buster

RUN pip install wheel
RUN pip install supervisor

COPY rpc/ rpc/
COPY docker-entrypoint.sh bin/

RUN cd rpc/ && python setup.py bdist_wheel && pip install dist/*

WORKDIR central/
COPY jaws_central/ jaws_central/
COPY requirements.txt .
COPY setup.py .
COPY MANIFEST.in .
COPY tests/ tests/
RUN python setup.py bdist_wheel && pip install dist/*

ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["supervisord -n"]
