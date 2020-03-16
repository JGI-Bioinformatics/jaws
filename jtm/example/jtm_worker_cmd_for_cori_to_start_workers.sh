source ~/venv/bin/activate && jtm-worker -wt dynamic -tp cori_worker_test -cl cori -c 32 -t 00:10:00 -m 1GB -C haswell -N 1 -nw 4 --qos genepool --account fungalp 
