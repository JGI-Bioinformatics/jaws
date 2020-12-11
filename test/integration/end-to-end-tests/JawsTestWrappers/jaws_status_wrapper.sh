#!/bin/bash
set -eo pipefail
# set -u: affects variables. When set, a reference to any variable you haven't previously defined 
# set -e: forces an exit with a return code of >0
# set -o pipefail: This setting prevents errors in a pipeline from being masked. If any command in a pipeline fails, that return code will be used as the return code of the whole pipeline. By default, the pipeline's return code is that of the last command - even if it succeeds.

ANALYSIS_FILE=
TABLE_NAME=
SITE=
JAWS_RELEASE=
RQC_VENV=/global/dna/projectdirs/PI/rqc/prod/versions/jgi-rqc-autoqc/config/rqc38.sh
AUTOQC_TOOL=/global/dna/projectdirs/PI/rqc/prod/versions/jgi-rqc-autoqc/bin/autoqc_tool_gp.py
JSON=/global/cscratch1/sd/jfroula/JAWS/jaws/test/integration/end-to-end-tests/TestWDLs/test.json
WDL=/global/cscratch1/sd/jfroula/JAWS/jaws/test/integration/end-to-end-tests/TestWDLs/test.wdl

usage()
{
cat<<EOF
  Usage: $0 
    <-r jaws release [jaws-prod|jaws-staging|jaws-dev] (required)> 
    <-a output analysis file name (required)> 
    <-t name of threashold file as defined in meta{} (required)>
    <-s compute site [cori|jgi] (required)>
EOF
exit 1
}


while getopts 'a:t:s:r:h' OPTION
do 
  case $OPTION in 
  a)    ANALYSIS_FILE="$OPTARG"
        ;;
  t)    TABLE_NAME="$OPTARG"
        ;;
  s)    SITE="$OPTARG"
        ;;
  r)    JAWS_RELEASE="$OPTARG"
        ;;
  h)    usage
        ;;
  ?)    usage
        ;;
  esac
done

# check for all required args
ARGS=$*
if [[ $ARGS =~ "-a" ]] && [[ $ARGS =~ "-t" ]] && [[ $ARGS =~ "-s" ]] && [[ $ARGS =~ "-r" ]]; then
    echo;
else
    echo -e "Missing some required arguments\n"
    usage
fi

echo $JAWS_RELEASE
if [[ $JAWS_RELEASE == "jaws-prod" ]] || [[ $JAWS_RELEASE == "jaws-staging" ]] || [[ $JAWS_RELEASE == "jaws-dev" ]]; then
	echo 
else
	echo "Please use one of the following for the site: [jaws-prod|jaws-staging|jaws-dev]."
fi

#TODO jeff. this should be removed once Stephan changes jaws to jaws-prod for "--profile" argument
if [[ $JAWS_RELEASE == "jaws-prod" ]]; then
	JAWS_RELEASE="jaws"
fi

# 
# Create the Analysis File
# 
source ~jfroula/jaws-prod.sh
echo "jaws_status_test.py -w $WDL -i $JSON -a $ANALYSIS_FILE -s $SITE -n $TABLE_NAME"
	  jaws_status_test.py -w $WDL -i $JSON -a $ANALYSIS_FILE -s $SITE -n $TABLE_NAME
deactivate

#
# Compare Analysis and Threshold files and upload results to GUI
#
source $RQC_VENV
echo "$AUTOQC_TOOL --profile $JAWS_RELEASE qc $ANALYSIS_FILE $TABLE_NAME"
      $AUTOQC_TOOL --profile $JAWS_RELEASE qc $ANALYSIS_FILE $TABLE_NAME
