#!/bin/bash
cd test/integration/end-to-end-tests/nightly-tests/health-robot
pytest -n 4 --capture=no --verbose test_fq_count_all_sites.py
