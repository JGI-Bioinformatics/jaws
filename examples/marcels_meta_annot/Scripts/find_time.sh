#!/bin/bash
i=$1
start=$(find call-s_annotate/shard-${i}/sa.s_annotate/*/call-trnascan/trnascan.trnascan/*/call-trnascan_ba/execution/ -name script -exec ls -l {} \; | awk '{print $6,$7,$8}'| head -1)
end=$(find call-f_annotate/shard-${i}/*/*/call-product_name/ -name script -exec ls -l {} \; | awk '{print $6,$7,$8}' | head -1)
echo "shard $i [$start : $end]"

