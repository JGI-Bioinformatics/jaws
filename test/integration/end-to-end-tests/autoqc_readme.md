# Summary
This repo should contain all the code for generating the analysis and threshold files to satisfy the AutoQC JAWS "scorecard" tests.

Here are some useful links:

* [Description of AutoQC](https://code.jgi.doe.gov/qaqc/autoqc/-/wikis/AutoQC-User-Help)
* [Example Threshold files](https://code.jgi.doe.gov/qaqc/autoqc-data)
* [Main repo](https://code.jgi.doe.gov/qaqc/autoqc) 
    see under thresholdValidator/ there will be the bin/auto-qc script and the auto-qc library dir

## Example commands for running tests through AutoQC 
Below are example steps to complete "scorecard" tests.  The process starts with a submission to JAWS and ends with a AutoQC report at https://rqc-10.jgi-psf.org/autoqc_gp/report/.  Two files are required, an Analysis file (yaml) and a Thresholds file (json (or yaml)).

The script to create the Analysis file has been adapted from Angie Kollmer's code "task_log_check_loop.py" and "submit_loop.py" that will submit a WDL to JAWS and wait for it to complete.  Then submit another test command, the results of which will be converted to an Analysis file. 

## Description of JawsTestWrappers scripts
#### jaws_status_test.py
```
This script essentially takes the raw json output from the JAWS cmd and converts the values to something   
that can be validated by the AutoQC evaluators:

  1) equals
  2) greater_than
  3) less_than
  4) greater_equal_than
  5) less_equal_than
  
  See [Description of AutoQC](https://code.jgi.doe.gov/qaqc/autoqc/-/wikis/AutoQC-User-Help) for more.

There are two tests functions in jaws_status_test.py now, as examples, 

   * test_upload_task_id
   * test_status

The functions convert the values of the raw JAWS output to a 1 or 0, for pass or fail.
```

#### parsing_functions.py
This is a module adapted from submit_loop.py.  It includeds functions to submit a jaws run and wait for it to complete.  You can submit one or more runs depending on the function you use.  It also has a function to create the Analysis yaml file.


## Run an example jaws submission
....first, clone this repo and cd into it.

Our example test will be to see that a JAWS submission was successfull.

## Create the Analysis file

```
. ~/jaws-prod.sh
./JawsTestWrappers/jaws_status_test.py -w TestWDLs/test.wdl -i TestWDLs/test.json -a analysis.yaml \
  -s cori -n jaws_submit_1
deactivate
```

You should now have an Analysis yaml file called `analysis.yaml` that looks like this:
```
metadata:
    name: jaws_submit_1
data:
    upload_task_id: 1
    status: 1
```




## Register the Threshold file

make sure you deactivate any previous environments with `deactivate` and then activate this venv:

```
.  /global/dna/projectdirs/PI/rqc/prod/versions/jgi-rqc-autoqc/config/rqc38.sh
```

Register the Threshold

```
/global/dna/projectdirs/PI/rqc/prod/versions/jgi-rqc-autoqc/bin/autoqc_tool_gp.py \
-p jaws add threshold ThresholdFiles/jaws-submit-threshold.json [--update]

# use the --update flag if you change the Threshold file but don't want to change the version.

# currently the threshold file can be reqistered as above, but an additional, manual step
# must be completed by Stephan before analysis files can be submitted against that threshold file
# (the manual step is creating a new collection in the mongo db)  
```

## Compare Analysis and Threshold files and upload results to GUI
```
/global/dna/projectdirs/PI/rqc/prod/versions/jgi-rqc-autoqc/bin/autoqc_tool_gp.py \
-p jaws qc analysis.yaml jaws_status

where:
  analysis.yaml = analysis file
  jaws_status = threshold name (defined in the threshold file metadata)
```

## See the AutoQC web table

```
https://rqc.jgi-psf.org/autoqc_gp/report/
click on `Jaws_status` and then the "Search" button.
```
