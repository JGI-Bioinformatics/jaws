## CHANGELOG

### Feb 2, 2021
* I added a wait_for_wdl function to the test functions in test_jaws_cmds.py. Now each function is reponsible for waiting for the WDL run to finish. This means I can use just the one fixture "submit_wdl" for both cases, for testing while a run is running, and tests for when the run is complete.
* the env" variable used for the fixtures is hard coded still (i.e. prod). We need to figure out how to make it a param. Passing it in as a param from the command line (i.e. pytest --env prod) doesn't work when the fixture scope="module".  If I remove the scope param to the fixtures, env is a passable variable, but each fixture because scoped to "function", meaning submit_wdl runs for each function; and thus we have too many WDL submissions.
* added task_run_stats.py which checks that JAWS runs are returning correct logging info from the cli
* added inputs_outputs_available.py which checks that the input WLD & json files are saved to the output dir, and that the raw cromwell dir is copied to the output dir.

### Jan 27, 2021
* first commit 
* includes jeff's health-check-robot code
* includes score-card-tests/test_jaws_cmds.py. See the README.md to run this pytest example.
