## Summary
This is an example running a pipeline in JAWS. 
It demonstrates the usage of a docker image, a subworkflow and how sharding of an input file is done.

## Set up JAWS Environment
```
# make sure you have the "activate" command
module load python  

# create this directory if you don't already have it
mkdir ~/.conda/envs

# create a local conda environment for jaws to make life easier in the future
ln -s /global/project/projectdirs/jaws/prod/cli/ ~/.conda/envs/jaws
source activate jaws 
... or use this instead if you have conda set up this way (you'll know if you do).
conda activate jaws

# verify you have the command
jaws
```

## Log in (do this once)
```
jaws login <NERSC password + MFA>
```

## Running
```
# clone the example code
git clone https://gitlab.com/jfroula/jaws-example-wdl.git

cd jaws-example-wdl/jaws-alignment-example

# run jaws submit <workflow> <inputs> <subworkflow>
jaws submit alignment.wdl inputs.align.json

# you should see something like this
Successfully queued job c6d403a0-bc07-495f-8250-7132593e7d7d
```

## monitoring the job
```
# make sure you copy the id of the job submission, if you didn't you can run this to see your run's id
jaws queue

# check jaws status
jaws status c6d403a0-bc07-495f-8250-7132593e7d7d
```

## output
check jaws results directory.
The last line of the json output `workflowRoot` should be a directory where the results are. 
You can copy and flatten the directory structure by using the `wfcopy.py` script that is available with the jaws environment.

```
# look for workflowRoot= near bottom of metadata output
jaws metadata c6d403a0-bc07-495f-8250-7132593e7d7d

# copy & simplify directory structure (including logs) in a results directory of your choosing.
wfcopy.py <workflowRoot> <your_results_dir>
```
