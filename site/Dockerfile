ARG JAWS_VERSION
FROM rpc:${JAWS_VERSION} as builder

ENV JAWS_VERSION=${JAWS_VERSION}
WORKDIR /usr/app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .

FROM builder as dev
RUN mkhomedir_helper jaws
RUN pip install -e .
USER jaws

FROM builder as prod
RUN python setup.py bdist_wheel
RUN pip install dist/*
USER jaws
