#!/bin/bash
INITIAL=$1
FINAL=$2
WDL=$3
if [[ ! $INITIAL ]] || [[ ! $FINAL ]] || [[ ! $WDL ]]; then
    echo "Usage: $0 <'initial\/image'> <'final\/image'> <wdl>"
    exit 1
fi

#sed 's/bfoster1\/img-omcis:0.0.7/jfroula\/img-omics:0.0.8/' genemark.wdl
#echo sed \"s/$INITIAL/$FINAL/\" $WDL
sed -i "s/$INITIAL/$FINAL/" $WDL
