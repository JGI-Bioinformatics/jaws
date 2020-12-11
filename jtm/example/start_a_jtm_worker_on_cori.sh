# Remove "--dry_run" for testing
jtm worker \
  -wt dynamic \
  -cl cori \
  -t 00:10:00 \
  -N 1 \
  -m 10G \
  -c 1 \
  -A fungalp \
  -q genepool_special \
  -p test_worker_cori \
  --dry_run
