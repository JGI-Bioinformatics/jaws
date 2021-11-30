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

