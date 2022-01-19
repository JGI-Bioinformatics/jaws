FROM python:3.8-buster

RUN pip install wheel

WORKDIR jaws
COPY jaws .
RUN cd rpc && python setup.py bdist_wheel && pip install dist/*

WORKDIR /
