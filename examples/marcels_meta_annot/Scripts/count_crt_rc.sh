#!/bin/bash
id=$1
cat /global/cscratch1/sd/jaws_jtm/prod/cromwell-executions/annotation/${id}/call-s_annotate/shard-*/sa.s_annotate/*/call-crt/crt.crt/*/call-run/execution/rc | wc -l
cat /global/cscratch1/sd/jaws_jtm/prod/cromwell-executions/annotation/${id}/call-s_annotate/shard-*/sa.s_annotate/*/call-crt/crt.crt/*/call-run/execution/rc | sort -u
