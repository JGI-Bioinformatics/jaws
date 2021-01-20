# End to End Integration Testing
This repo should contain all the code for generating the analysis and threshold files to satisfy the AutoQC JAWS "scorecard" tests.

Here are some useful links:

* [Description of AutoQC](https://code.jgi.doe.gov/qaqc/autoqc/-/wikis/AutoQC-User-Help)

* [Example Threshold files](https://code.jgi.doe.gov/qaqc/autoqc-data)

* [Main repo](https://code.jgi.doe.gov/qaqc/autoqc)
   see under thresholdValidator/ there will be the bin/auto-qc script and the auto-qc library dir

## Example commands for running tests through AutoQC
Below are the steps to complete the "scorecard" tests.  The tests can be viewed as an AutoQC report at https://rqc.jgi-psf.org/autoqc_gp/report/.  You need to generate two files, 
1) an AutoQC report file (yaml) and 
2) a thresholds file (json).

There are two template files that you can use:

* JawsTestWrappers/template.py file that should be copied and used for your tests and will create the report yaml file.
* ThresholdFiles/template.json can be used to generate a corresponding threasholds file.

## Procedure
1) What I do is first write the tests (e.g. JawsTestWrappers/test_jaws_cmds.py). This will produce an autoQC file that is used for the GUI, but also creates an intermediate file, e.g. analysis.yaml.

   The tests can be run like so
   ```
   JawsTestWrappers/test_jaws_cmds.py -w TestsWDLs/fq_count.wdl -i TestsWDLs/fq_count.json -s cori -e prod
   ```
   A log file is also create (e.g. test_jaws_cmds.log).  
2) You can model the thresholds file from the analysis.yaml file, the path to which can be found in the log file.
3) The thresholds file needs to be registered into the AutoQC mongo database. This can be done by
   ```
   source /global/dna/projectdirs/PI/rqc/prod/jgi-rqc-autoqc/config/rqc38.sh

   # There are different profiles (jaws_prod, jaws_staging, and jaws_dev).
   # use the --update flag if you change the Threshold file but don't want to change the version.
   autoqc_tool_gp.py --profile jaws_prod add threshold test_jaws_cmds.json [--update]
   ```

4) view report https://rqc.jgi-psf.org/autoqc_gp/report/

## Designing the Thresholds File
These are AutoQC terms that you can use in the thresholds file.
```
1) equals
2) greater_than
3) less_than
4) greater_equal_than
5) less_equal_than
```
See [Description of AutoQC](https://code.jgi.doe.gov/qaqc/autoqc/-/wikis/AutoQC-User-Help) for more.

## Library Functions
#### parsing_functions.py

This is a module adapted from submit_loop.py.  It includeds functions to submit a jaws run and wait for it to complete. Please add shared functions here. 
