v1.1
- Change to manage resources in four types: small, medium, large, xlarge, mainly based on memory spec
- User per-site configuration 
- Request compute node based on the memory and cpu requests
- Add `stop_worker.sh` to the SLURM job script to forcefully terminate Condor worker processes so that `condor_status` doesn't show the slot(s) for the compute node
- Increase the default memory request size from 512MB to 10GB for all sites
- Simplify the condor central and worker configs and use password based authentication


 
v1.0 
- Request compute node based on the number of idle Condor jobs