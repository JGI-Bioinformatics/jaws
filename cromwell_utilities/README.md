# Cromwell Utilities

## Installation

    cd ~
    git clone https://code.jgi.doe.gov/advanced-analysis/jaws.git
    python3 -m venv ./cromutilenv
    source ~/cromutilenv/bin/activate
    cd ./jaws/cromwell_utilities
    python setup.py develop
    export CROMWELL_URL="localhost:50010"
    cromwell-utils --help
