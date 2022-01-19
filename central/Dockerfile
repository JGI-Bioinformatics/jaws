ARG JAWS_VERSION
FROM jaws_rpc:${JAWS_VERSION}

RUN pip install wheel

COPY docker-entrypoint.sh bin/

WORKDIR app
COPY jaws/central .
RUN python setup.py bdist_wheel && pip install dist/*

ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["jaws-central --log /tmp/central-rest.log --config /etc/jaws-central.conf --log-level DEBUG rest"]
