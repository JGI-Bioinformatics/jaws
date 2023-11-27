FROM python:3.10-slim

# Install security updates, and some useful packages.
#
# * Make sure apt-get doesn't run in interactive mode.
# * Update system packages.
# * Pre-install some useful tools.
# * Minimize system package installation.
RUN export DEBIAN_FRONTEND=noninteractive && \
  apt-get update && \
  apt-get -y upgrade && \
  apt-get install -y --no-install-recommends tini procps npm \
  build-essential git make zip && \
  apt-get -y clean && \
  rm -rf /var/lib/apt/lists/*

# Install requirements
WORKDIR /code

COPY ./requirements.txt /code/requirements.txt

RUN pip install --no-cache-dir -r /code/requirements.txt

# Add repository code
COPY . /code
RUN pip install --no-cache-dir --editable .
# Prepare for C crashes.
ENV PYTHONFAULTHANDLER=1

CMD ["npx" "exec" "nodemon" "--exec" "jaws-site" "-h"]
