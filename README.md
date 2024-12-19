# Documentation for JGI Analysis Workflow Service (JAWS)

[![pipeline status](https://code.jgi.doe.gov/advanced-analysis/jaws-site/badges/main/pipeline.svg)](https://code.jgi.doe.gov/advanced-analysis/jaws-site/commits/main) [![coverage report](https://code.jgi.doe.gov/advanced-analysis/jaws-site/badges/main/coverage.svg)](https://code.jgi.doe.gov/advanced-analysis/jaws-site/commits/main)

### Resources for JAWS Users
Full [documentation](https://jaws-docs.readthedocs.io) for running and installing JAWS is located here.


# JAWS Site Services
The JAWS Site Services is a Python-based system that facilitates the management 
of file transfers, monitors and logs Cromwell WDL (Workflow Description Language) workflows, and handles messaging 
between JAWS Central and compute sites. Designed for scalability and reliability, it ensures seamless communication 
and processing for distributed bioinformatics workflows.

JAWS Site is our name for a collection of services that consist of the following:
- RPC server (Handles RPC messages from JAWS Centra (Handles RPC messages from JAWS Centra (Handles RPC messages from JAWS Centra (Handles RPC messages from JAWS Central))))
- Run daemon (monitors Cromwell runs and their status)
- Transfer daemon (monitors and transfers data either locally or via gLOBUS)
- Pool manager daemon (monitors condor queue for tasks and submits to Slurm)
- Task logger message consumer (keeps track of task logs)

## Getting Started
Prerequisites: 
- Python 3.10 >=, 
- RabbitMQ
- MySQL
- [Cromwell](https://github.com/broadinstitute/cromwell)
- [PDM](https://pdm-project.org/latest/)

### Installation
1. Clone Repository 
    `git clone https://code.jgi.doe.gov/advanced-analysis/jaws-site.git`
2. Install dependencies
    `pdm install`  
    This will: 
      - Install all dependencies specified in the pyproject.toml file.
      - Create and use a virtual environment automatically (if PDM’s virtualenv support is enabled).
3. Verify installation
   `jaws-site -h`

   
## Contributing
Developers
* Edward Kirton
* Mario Melara
* Elais Player
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

## License
This project is licensed under the [BSD License](LICENSE)
