# Documentation for JGI Analysis Workflow Service (JAWS)

[![pipeline status](https://code.jgi.doe.gov/advanced-analysis/jaws/badges/dev/pipeline.svg)](https://code.jgi.doe.gov/advanced-analysis/jaws/commits/dev) [![coverage report](https://code.jgi.doe.gov/advanced-analysis/jaws/badges/dev/coverage.svg)](https://code.jgi.doe.gov/advanced-analysis/jaws/commits/dev)

### Resources for JAWS Users
Full [documentation](https://jaws-docs.readthedocs.io) for running and installing JAWS is located here.


### Resources for JAWS Developers

### Installing JAWS site
See docs/install_Cli_and_Site.md  
See docs/startingServices.md  


### Requirements

An AMQP service (e.g. RabbitMQ) and relational database (e.g. MySQL) are required.  Central and each Site may have their own.  Services don't cross-talk to other databases.

#### RabbitMQ

Requires one user "jaws" that is used for all deployments (e.g. dev/staging/prod).  Deployments are given separate namespaces by defining distinct vhosts:

  - jaws_(DEPLOYMENT_NAME)
  - jtm_(DEPLOYMENT_NAME)

e.g jaws_dev, jtm_dev, jaws_staging, jtm_staging, jaws_prod, jtm_prod

Queues are automatically defined by the deployment scripts and don't need to be created in the RMQ admin pages.  Only the user and vhosts needs to be created outside of the JAWS CI/CD.


#### MySQL

Requires one user "jaws" that is used for all deployments (e.g. dev/staging/prod).  There are separate databases for each service/deployment:

  - cromwell_(SITE_NAME)_(DEPLOYMENT_NAME)
  - jaws_(SITE_NAME)_(DEPLOYMENT_NAME)
  - jtm_(SITE_NAME)_(DEPLOYMENT_NAME)

e.g. jaws_cori_dev, cromwell_cori_dev, jtm_cori_dev, ...

Additional, jaws-central requires:

  - jaws_central_(DEPLOYMENT_NAME)

e.g. jaws_central_dev, jaws_central_staging, jaws_central_prod


### Installing all services 
You can install using two modes, developer mode and build mode:
you need to know when to you develop vs build.  
	develop => instant changes since it's just a symlink to your src
	build => the installed code will not change if src changes. This mode builds everything so that you will able to install the source using pip into your venv

if you're editing the code (i.e. developing) and wish to test your changes, it's inconvenient to reinstall via pip for every change, so you use develop.  when done and you're satisfied with the changes, use build to install a copy of your package.

For an install on your mac, everything is installed via develop.  So if you switch branches, all the servers will be running whatever version of the software is in that branch.  for shared installations, you'd want to use develop only while you're tweaking things, install via build before you leave because someone else may change the branch in the src dir.

Example commands using develop:

```
git clone https://code.jgi.doe.gov/advanced-analysis/jaws.git
module load python/3.8-anaconda-2020.11
python3 -m venv ~/venv/jaws-test
source ~/venv/jaws-test/bin/activate
pip install black
pip install flake8
pip install pytest
pip install pytest-cov
pip install pytz

# inside the jaws repo 
# RPC (do this first)
pip install -r rpc/requirements.txt 
cd rpc && python setup.py develop

# Client
cd ../ && pip install -r client/requirements.txt
cd client && python setup.py develop
export JAWS_CLIENT_CONFIG=/global/cfs/projectdirs/jaws/jaws-dev/jaws-dev.conf

# Central
cd ../ && pip install -r central/requirements.txt
cd central && python setup.py develop

# Site
cd ../ && pip install -r site/requirements.txt
cd site && python setup.py develop

# JTM
cd ../ && pip install -r jtm/requirements.txt
cd jtm && python setup.py develop

# create a file called womtool in ~/venv/jaws-test/bin
 
  #!/bin/bash
  java -jar /global/cfs/projectdirs/jaws/cromwell/womtool.jar $*

chmod 755 ~/venv/jaws-test/bin/womtool

# Add your credentials. The jaws.conf file will contain your private jaws token
export JAWS_CLIENT_CONFIG=/global/cfs/projectdirs/jaws/jaws-dev/jaws-dev.conf
export JAWS_USER_CONFIG=~/jaws.conf

# test that environment is set up correctly
make test

# To deactivate the venv:
deactivate
```

Example commands using build:

```
git clone https://code.jgi.doe.gov/advanced-analysis/jaws.git
# repeat the same steps as above, except use "build" instead of "develop"
#   python setup.py develop

```

### cromwell-utils
Example installing cromwell-utils using build mode

```
# do this once
module load python/3.7
python -m venv ~/cromvenv
cd <path_to_repo>/jaws/cromwell_utilities/
python setup.py [develop|build]
pip install .  # run this if you used build in last step

# now do this every time
export TMPDIR=/"global/scratch/$USER"
export CROMWELL_URL=localhost:50010
source ~/cromvenv/bin/activate
cromwell-utils
```

### docker-compose services 
There is a `docker-compose.yml` setup for local deployments. To run the services be sure to have the db, rabbitmq services up first:
```console
docker-compose up -d central_db
docker-compose up -d central_rabbitmq
docker-compose up -d site_db
docker-compose up -d site_rabbitmq
```  

This will bring up all the required external services. You will then want to bring up the rest of the services: 

```console
docker-compose up -d
```

#### local development cycle 
If you want to develop jaws locally, you will need to re-build the images to make sure your current code is uploaded. 
The Dockerfiles copy the source code into the container, and then create all of the configuration scripts. Supervisord
is instaled in the containers and will start and run the services.  

There is a `./build.sh` script used to help create the docker images.

```console

$ ./build.sh
Usage: ./build.sh <jaws_service> <version>
```

The version tag can be customized and is used to simply tag your version of the docker images. Note that if you update
the version tag you will need to also update the docker-compose.yml. 

Example:  

```console
$ ./build.sh central 2.0
Sending build context to Docker daemon  189.4kB
Step 1/15 : FROM python:3.8-buster
 ---> 28a4c88cdbbf
Step 2/15 : RUN pip install wheel
 ---> Using cache
 ---> 7528d629c501
Step 3/15 : RUN pip install supervisor
 ---> Using cache
 ---> 39fb4b3d2e65
Step 4/15 : COPY rpc/ rpc/
 ---> Using cache
 ---> a814ea36fe3e
Step 5/15 : COPY docker-entrypoint.sh bin/
 ---> Using cache
 ---> 35828b47773d
Step 6/15 : RUN cd rpc/ && python setup.py bdist_wheel && pip install dist/*
 ---> Using cache
 ---> 9a62f1e6ca82
Step 7/15 : WORKDIR central/
 ---> Using cache
 ---> ceb7a9f4a08a
Step 8/15 : COPY jaws_central/ jaws_central/
 ---> Using cache
 ---> a939e24c8572
Step 9/15 : COPY requirements.txt .
 ---> Using cache
 ---> 8d5f1cae6d94
Step 10/15 : COPY setup.py .
 ---> Using cache
 ---> 0632054b56bc
Step 11/15 : COPY MANIFEST.in .
 ---> Using cache
 ---> 6f4e0cfc2d8a
Step 12/15 : COPY tests/ tests/
 ---> Using cache
 ---> 38ead4b5cb8b
Step 13/15 : RUN python setup.py bdist_wheel && pip install dist/*
 ---> Using cache
 ---> 354b34be6cbb
Step 14/15 : ENTRYPOINT ["docker-entrypoint.sh"]
 ---> Using cache
 ---> 6daae6b070a2
Step 15/15 : CMD ["supervisord -n"]
 ---> Using cache
 ---> 808fa3c83056
Successfully built 808fa3c83056
Successfully tagged jaws_central:2.0
```  

There is a Dockerfile contained within each subdirectory along with a docker-entrypoint.
The entrypoint scripts are meant to create a configuration file and populate the fields 
from environment variables. Each container uses supervisord to manage the multiple processes
needed for each application.  

## Contributing
Developers
* Edward Kirton
* Mario Melara
* Georg Rath
* Seung-jin Sul
* Stephan Trong
* Kelly Rowland
* Nick Tyler

Documentation and WDL authoring
* Jeff Froula
 
User support
* Jeff Froula
* Ed Kirton
* Ramani Kothadia
* Seung-Jin Sul

System and integration testing
* Jeff Froula
* Angie Kollmer
* Ramani Kothadia

Functional design consulting, project and resource coordination
* Kjiersten Fagnan
* Hugh Salamon

