FROM python:3.8-buster

RUN pip install wheel
RUN pip install supervisor
COPY docker-entrypoint.sh bin/

COPY rpc/ rpc/
RUN cd rpc/ && python setup.py bdist_wheel && pip install dist/*

WORKDIR site/
COPY jaws_site/ jaws_site/
COPY requirements.txt .
COPY setup.py .
COPY sql_tables sql_tables/
COPY tests/ tests/
RUN python setup.py bdist_wheel && pip install dist/*

ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["supervisord -n"]
