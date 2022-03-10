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

## Local Development
To deploy the JAWS services (except for the backend) locally you will want to make use
of the docker-compose.yml file and modify as needed. Some volume mounts are required in order for file uploads to
be processed by the containers. You will want to set your `DATA_HOME` with the locations of your
configuration files, log files and data upload files. Here is an example: 


```console
❯ pwd
/Users/mamelara/data/jaws
~/data/jaws on ☁️  (us-west-2)
❯ tree -L 2
.
├── configs
│   ├── auth
│   ├── central
│   ├── central-rpc
│   ├── cromwell
│   ├── daemon
│   └── site-central
├── laptop
│   ├── hello_world
│   └── output
├── local
│   ├── cromwell-workflow-logs
│   └── uploads
└── logs
    ├── auth
    ├── central
    ├── central-rpc
    ├── cromwell
    ├── cromwell-workflow-logs
    ├── daemon
    ├── jaws-auth
    └── site-central

22 directories, 0 files
> export DATA_HOME=/Users/mamelara/data/jaws
```

The important thing to note is the name of the files, though these can be easily changed
in the `docker-compose.yml` file if you wish to use a different naming scheme. One you've
created the following directory tree, the next step is to build the images.

The `laptop` and `local` directories are just there to "mock" locations on different filesystems. You can name these whatever directory
you want as long as you make sure the volume mounts match the names.

### Image building
At the root of the JAWS directory is a `build.sh` script that takes arguments for 
 the service name, the version and environment (eg - dev, prod). If you want to build central you will
want to run the following command: 

```console
> ./build.sh central 2.6.0 dev
```

This will build the development environment for your container script. This is useful if you want the code you are actively working on to be deployed in the container
Because rpc is a dependency for central, and site services, you will first want to build rpc docker container. 

```console
> ./build.sh rpc 2.6.0 dev
```

Once your images are built you can then modify the `docker-compose.yml` with the apppropriate 
image tags and then deploy the services. You'll want to deploy the mysql and rabbitmq services
first: 

```console
> docker-compose up -d db rabbitmq
```
You will want to setup your database first by creating some databases in mysql: 
You can enter the container either with `docker exec` or with the docker desktop. It is recommended 
to use the docker desktop since it provides a nice GUI. If you want to
```console
> docker exec -it <CONTAINER ID or CONTAINER NAMES> /bin/bash
```

If you're not sure of the container id or the container name you can run `docker ps`.
You will want to create the databases inside mysql: 

```console
$ mysql -u root -p 
$ Password: ********
$ CREATE DATABASE IF NOT EXISTS jaws_cromwell_dev;
$ CREATE DATABASE IF NOT EXISTS jaws_local_dev;
$ GRANT ALL PRIVILEGES ON * . * TO 'jaws'@'%';
```

Next you can bring up the other services without specifying their names.

```console
> docker-compose up -d
```

You should now have a deployed version of JAWS (except for the backend) on your local dev environment. To confirm
visit: http://localhost:80/api/v2/status and you should see the JAWS service statuses displayed in your browser. 

In order to use this local instance you will want to create a user in the mysql container. Enter the container and then
open up the mysql console. There is a table in the `jaws_central_dev` database called users. You will want to enter your
user information there. 

````console
$ mysql -u jaws -p
$ Password:
$ use jaws_central_dev;
$ insert into users (id, email, name, is_admin, is_dashboard, jaws_token) values ("yourusername", "youremail", true, false, "yourtoken");1
````

This should insert the first user in the database which will be your user. 

### Local Development Flow
The real power of deploying locally comes in being able to have your own instance of JAWS available. You have
all other services available for you with little setup. Let's look at an example on how to develop JAWS central code 
locally.

You'll want to use the `./build.sh` script to make sure you build the dev container version.

```console
./build.sh central 2.6.0 dev
```

This will build the central dev image with the tag `2.6.0-dev`. Next you will want to make sure you bind mount the source code you are working on: 

```yaml

volumes:
	- $DATA_HOME/configs/central/jaws-central.conf:/etc/config/jaws-central.conf
	- $DATA_HOME/laptop:/data/jaws/laptop:rw
	- $DATA_HOME/local:/data/jaws/local:rw
	- $DATA_HOME/logs/central:/var/log/jaws-central:rw
	- ./central:/usr/app
```

Note that the working directory is located in `/usr/app`. This is where the source code will be installed in the container.

You'll also want to just run an interactive shell without any commands running:

```yaml
central:
	image: central:2.6.0-dev
	stdin_open: true  # add these to the docker-compose file
	tty: true
```

Then you can run a `docker-compose up -d`. And your container will be running. You can enter the shell once again
with `docker exec -it <CONTAINER_ID> /bin/bash`. Once you are in there you can install your package inside
the container by running `pip install -e .`

Since you mounted your code, any changes reflected in your code will be reflected in the container. You can then
run that service with all the updated code changes you made to do some live testing. 

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
* Ramani Kothadia
 
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
* Stephan Chan

