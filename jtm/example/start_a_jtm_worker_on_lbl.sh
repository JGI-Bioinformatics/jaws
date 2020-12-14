# Remove "--dry_run" for testing
jtm worker \
  -wt dynamic \
  -cl lbl \
  -t 00:10:00 \
  -N 1 \
  -m 10G \
  -c 1 \
  -A lr_jgicloud \
  -q condo_jgicloud \
  -p test_worker_lbl \
  --dry_run
