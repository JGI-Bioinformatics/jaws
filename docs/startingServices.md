# Start JAWS Services

## Installing Services for JAWS 
How to Start up the Required Services for JAWS

## JAWS Client
```
ssh cori.nersc.gov
source /global/cfs/projectdirs/jaws/miniconda3/bin/activate /global/cfs/projectdirs/jaws/cli/jaws-dev/
jaws util site
```

## JAWS Central
```
ssh jaws.lbl.gov   (no MFA)
source /opt/miniconda3/bin/activate /opt/jaws/dev/
cd /opt/jaws/dev/opt/central
./server.py &
````
NOTE: this will likely be renamed to jaws-central

## JAWS Auth
```
ssh jaws.lbl.gov   (no MFA)
source /opt/miniconda3/bin/activate /opt/jaws/dev/
cd /opt/jaws/dev/opt/auth
./server.py &
```
NOTE: this will likely be renamed to jaws-auth

## JAWS Site-NERSC
```
ssh cori.nersc.gov
ssh cori20
collabsu jaws   (no MFA)
source /global/cfs/projectdirs/jaws/miniconda3/bin/activate /global/cfs/projectdirs/jaws/site/dev/
server.py
```
NOTE: this will likely be renamed to jaws-site

## JAWS Site-LBNL (NEW)
```
ssh lrc-services.lbl.gov
sudo systemctl start jawsdev.service
```

## JAWS Site-LBNL (OLD)
```
ssh lrc-services.lbl.gov
source /global/home/groups-sw/lr_jgicloud/miniconda3/bin/activate /global/home/groups-sw/lr_jgicloud/site/dev/
server.py
```
NOTE: this will likely be renamed to jaws-site

# Cromwell-NERSC
```
ssh cori.nersc.gov
ssh cori20
collabsu jaws_jtm   (no MFA)
source /global/cfs/projectdirs/jaws_jtm/anaconda3/bin/activate /global/cfs/projectdirs/jaws_jtm/dev/jtm
start_cromwell.sh
```

## Cromwell-LBNL
```
ssh lrc-services.lbl.gov
source /global/home/groups-sw/lr_jgicloud/anaconda3/bin/activate /global/home/groups-sw/lr_jgicloud/dev/jtm
start_cromwell.sh
```

## MySQL
```
ssh lrc-services.lbl.gov
ssh jaws-db.lbl.gov   (no MFA)
sudo systemctl start mysql.service
````

JAWS-Dashboard
TODO

--------------------------------------
## JTM-Manager-NERSC
#### On cori20 with jaws
```
collabsu jaws

# DEV
source /global/cfs/projectdirs/jaws/jtm/venv-dev/bin/activate && jgi-task-manager

# PROD 
source /global/cfs/projectdirs/jaws/jtm/venv/bin/activate && jgi-task-manager
```

####  On cori20 with jaws_jtm
```
collabsu jaws_jtm

# DEV
source /global/cfs/projectdirs/jaws_jtm/anaconda3/bin/activate /global/cfs/projectdirs/jaws_jtm/dev/jtm && jgi-task-manager

# PROD
source /global/cfs/projectdirs/jaws_jtm/anaconda3/bin/activate /global/cfs/projectdirs/jaws_jtm/prod/jtm && jgi-task-manager
```

Note: workers can be started with the same env activation + 
          `nohup jtm-worker -tp <user_pool_name> &`

## JTM-Manager-LBNL
By systemctl
```
sudo systemctl start/stop/status jtmdev.service
sudo systemctl start/stop/status jtmprod.service
```

Starting as user (i.e. as you).
```
# DEV
source /global/home/groups-sw/lr_jgicloud/anaconda3/bin/activate /global/home/groups-sw/lr_jgicloud/dev/jtm && jgi-task-manager

# PROD
source /global/home/groups-sw/lr_jgicloud/anaconda3/bin/activate /global/home/groups-sw/lr_jgicloud/prod/jtm && jgi-task-manager
```

Note: workers can be started with the same env activation + 
          `nohup jtm-worker -tp <user_pool_name> &`


