# Documentation for JGI Analysis Workflow Service (JAWS)

[![pipeline status](https://code.jgi.doe.gov/advanced-analysis/jaws-site/badges/main/pipeline.svg)](https://code.jgi.doe.gov/advanced-analysis/jaws-site/commits/main) [![coverage report](https://code.jgi.doe.gov/advanced-analysis/jaws-site/badges/main/coverage.svg)](https://code.jgi.doe.gov/advanced-analysis/jaws-site/commits/main)

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

Please see the detailed instructions in the local_development.md file in this repository.

## Pip-tools and pyproject.toml
`pyproject.toml` is the latest standard in configuring Python projects. JAWS has moved over to this over
`setup.py` but installing this project remains the same. There are two requirements files, `dev-requirements.txt` and
`requirements.txt`, which contain the pinned versions of JAWS dependencies. These two files were generated using [pip-tools](https://github.com/jazzband/pip-tools) as recommended by the [Python Packaging Guide](https://packaging.python.org/en/latest/guides/tool-recommendations/).

To generate these files you will want to first install pip-tools and then run `pip-compile -o requirements.txt pyproject.toml`.
The `pip-compile` command will then resolve the dependencies declared in `pyproject.toml` into the output file. You can
install these dependencies either using `pip install -r requirements.txt` or using `pip-sync`

Since setup.py is removed, you can no longer run `python setup.py install`. These calls have been deprecated in favor
of installing projects with `pip install .` There is an interesting discussion in this [post](https://blog.ganssle.io/articles/2021/10/setup-py-deprecated.html).


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

