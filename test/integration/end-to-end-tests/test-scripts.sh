#!/bin/bash

HOMEDIR=/global/cscratch1/sd/jfroula/JAWS/jaws/examples
REFDIR=referencing_db_and_shifter
SUBDIR=jaws-alignment-example

set -x
### testcase1
#cd $HOMEDIR/referencing_db_and_shifter
#echo "running test in $(pwd)"
#
#echo -n "testcase1..."
#jaws run submit $HOMEDIR/$REFDIR/test.wdl $HOMEDIR/$REFDIR/inputs.json $(pwd)/out nersc > out1 2>&1
#id=$(grep '"run_id":' out1 | awk -F':' '{print $2}' | tr -d ' ,')
#if [[ $id > 0 ]]; then
#	echo "Success: found id $id"
#else
#	>&2 echo "Failed: couldnt find id"
#fi

id=131
### testcase2
echo "testcase2..."
jaws run status $id > status1 2>&1
jaws run tasks $id > tasks1 2>&1
jaws run queue > queue1 2>&1

# test for state
echo -n "    state..."
state=$(grep status status1 | tr -d '" ,"' | awk -F':' '{print $2}')
list_of_possibles="uploading submitted running succeeded ready downloading finished"
if echo $list_of_possibles | grep -w $state > /dev/null; then
	echo "success: state is $state"
else
	>&2 echo "failed"
fi

# test for update time
echo -n "    updated..."
updated=$(grep updated status1 | tr -d '" ,"' | awk -F':' '{print $2}')
OLDIFS=$IFS
IFS='-' read -ra TIME <<< "$updated" # str is read into an array as tokens separated by IFS
if [[ ${#TIME[@]} == 3 ]]; then
	echo "success: updated $updated"
else
	>&2 echo "failed"
fi

# test for when entered & exited tasks
#[
#    [
#        "runblastplus_sub.task1",
#        "Done",
#        "2020-05-06T20:53:05.108Z",
#        "2020-05-06T20:54:28.898Z"
#    ],
#    [
#        "runblastplus_sub.task2",
#        "Done",
#        "2020-05-06T20:54:30.787Z",
#        "2020-05-06T20:54:32.894Z"
#    ]
#]

echo -n "    task interval..."
time1=$(head -5 tasks1 | tail -1 | tr -d '", ')
time2=$(head -6 tasks1 | tail -1 | tr -d '", ')
IFS='-' read -ra TIME1 <<< "$time1"
IFS='-' read -ra TIME2 <<< "$time2"
echo mytime ${TIME2[@]}
if [[ ${#TIME1[@]} == 3 ]] && [[ ${#TIME2[@]} == 3 ]]; then
	echo "success: found correct format for start time and end time"
else
	>&2 echo "failed"
fi

# A list of tasks with unique IDs for the run
echo -n "    number of task ..."
tasks=$(cat tasks1 | grep runblastplus_sub | tr -d '", ' | wc -l)
if [[ $tasks == 3 ]]; then
	echo "success: found $tasks tasks" 
else
	>&2 echo "failed"
fi

# Input,Output & WDL data for the run (paths)
echo -n "    input_file"
input_file=$(tail -17 referencing_db_and_shifter/queue1 | tr -d '", ' | grep input_file | awk -F':' '{print $2}')
if [[ -e $input_file ]]; then
	echo "success: found input_file $(basename $input_file)"
else
	>&2 echo "failed"
fi

echo -n "    output_dir"
output_dir=$(tail -17 referencing_db_and_shifter/queue1 | tr -d '", ' | grep output_dir | awk -F':' '{print $2}')
if [[ $output_dir ]]; then
	echo "success: found output_dir $(basename $output_dir)"
else
	>&2 echo "failed"
fi

echo -n "    wdl_file"
wdl_file=$(tail -17 referencing_db_and_shifter/queue1 | tr -d '", ' | grep wdl_file | awk -F':' '{print $2}')
if [[ -e $wdl_file ]]; then
	echo "success: found wdl_file $(basename $wdl_file)"
else
	>&2 echo "failed"
fi

### testcase 3
#With a Task-ID the user should be able to receive the following information through the JAWS CLI about the task:
#  stderr/stdout of the task
#  State
#  Information since when task is in said state
#  A list of transitions between states (entered state at time, left state at time)
#  Shell script as executed by JAWS
jaws run errors $id > errors1 2>&1
#jaws run submit TestsWDLs/test.wdl inputs.json $(pwd)/out nersc


### testcase 15


### testcase 20
#shared: Test that runtime {shared: 1} works.
#The shared: 1 option allows jaws user to re-use pools when they run another workflow and pools are still up.  The shared: 0 option only allows tasks within a run instance to share the pool.
#One thing to note is that if the workers removed by scancel, jtm doesn’t know it’s removed.
#once remove-pool is added to jaws command, you can use the command to remove workers.
#
#TESTCASE-21
#constraint: what happens if you set constraint: "knl" and run on lbnl?  NERSC options are knl and haswell. LBNL has only haswell.
#
#TESTCASE-22
#mem: what if this is out of bounds for the possible site resources. NERSC has limit of 450G machines.
#
#TESTCASE-23
#time: what if this is out of bounds for the possible site resources. NERSC has limit of 72hrs.
#
#TESTCASE-24
#poolname: what if you use small,med,large,xlarge
#
#TESTCASE-25
#node: limit ?
#
#TESTCASE-26
#nwpn: limit should be less than number of threads?
#
#TESTCASE-27
#subworkflows
#
