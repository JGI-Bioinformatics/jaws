## CHANGELOG

### Feb 25, 2021 
Jeff Froula
* added more tests
	- test_error_handling.py
	- test_inputs_outputs_available.py
	- test_resources_avail.py

### Feb 17, 2021
Jeff Froula
* converted the following to the new system:
    inputs_outputs_available.py
    subworkflow_behavior.py
    task_run_stats.py
* created a template.py file to be used for new tests.
* moved all test wdls and input json files into a WDLs folder.

### Feb 16, 2021
Jeff Froula
* moved the wait_for_wdl function into the submission_utils.py module so it can be shared code.

### Feb 12, 2021
Angie Kollmer
* found a way to have command line arguments like "env" and "site" and the use of fixtures. This allows us to only have a WDL submitted once per "module" or "session" and control by command line, the environment and which site the WDL is submitted to. 
* created a submission_utils.py file that have non-fixture functions that submit & wait for WDLs. These functions are called by the fixtures in conftest.py that have hard-coded values (i.e. wdl, input.json) but also capture the command line arguments (i.e. env & site).
* converted the test_jaws_cmds.py script to comply with her new system.

### Feb 9, 2021
Jeff Froula
* added subworkflow_behavior.py that checks JAWS handles subworkflows correctly.
* The conftest.py file now requires a ini file (i.e. fq_count.ini) to pass variables (ENV,WDL,INPUT_JSON,OUTDIR,SITE).  These variables are used in the submission of JAWS runs. The ini file is itself an environmental variable so when running pytest, you need to make sure MYINI_FILE is set.For example:  MYINI_FILE=subwdl.ini pytest --capture=no --verbose subworkflow_behavior.py

### Feb 2, 2021
* added a wait_for_wdl function to the test functions in test_jaws_cmds.py. Now each function is reponsible for waiting for the WDL run to finish. This means I can use just the one fixture "submit_wdl" for both cases, for testing while a run is running, and tests for when the run is complete.
* the env variable used for the fixtures is hard coded still (i.e. prod). We need to figure out how to make it a param. Passing it in as a param from the command line (i.e. pytest --env prod) does not work when the fixture scope="module".  If I remove the scope param to the fixtures, env is a passable variable, but each fixture because scoped to "function", meaning submit_wdl runs for each function; and thus we have too many WDL submissions.
* added task_run_stats.py which checks that JAWS runs are returning correct logging info from the cli
* added inputs_outputs_available.py which checks that the input WLD & json files are saved to the output dir, and that the raw cromwell dir is copied to the output dir.

### Jan 27, 2021
* first commit 
* includes jeffs health-check-robot code
* includes score-card-tests/test_jaws_cmds.py. See the README.md to run this pytest example.

