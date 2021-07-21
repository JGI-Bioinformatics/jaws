#!/bin/bash 
HERE=`pwd`
echo We are starting here: $HERE
cd test/integration/end-to-end-tests 
./single-test/test2.sh
