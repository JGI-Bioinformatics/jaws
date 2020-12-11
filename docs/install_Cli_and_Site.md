# JGI Analysis Workflow Service (JAWS)

## Installing JAWS site
JAWS Site (RPC server)
NOTE: you need a JAWS Site installation for every computing resource.  
It is an RPC (remote procedure call) Server which relays commands from the CLI to Cromwell via AMQP, as mediated by JAWS Central.

Alternate installation methods exist but Anaconda or miniconda is the easiest.
(JAWS Background)[https://jaws-docs.readthedocs.io/en/latest/Intro/how_jaws.html]  

## JAWS Requirements
Requirements:
It is assumed these services are already installed
* Globus
* MySQL (or other RDb server)
* RabbitMQ (or other AMQP server)
To install on your server:
      You need to be "jaws_jtm" user . 
      
```
collabsu jaws_jtm
export SITE=dev
```

Used only for steps below, repeat with “prod” for production branch
```
export SITE_BASE=/global/cfs/projectdirs/jaws_jtm/$SITE/site
```

Install (Anaconda)[https://www.anaconda.com/downloads] or miniconda3

```
conda create -p $SITE_BASE python=3
cd $SITE_BASE
mkdir -p opt/cromwell etc/conda/activate.d
echo export JAWS_SITE_CONFIG=$SITE_BASE/etc/jaws_client.ini  > etc/conda/activate.d/env_vars.sh
```

Get JAWS Site
```
cd $SITE_BASE/opt
git clone https://gitlab.com/jgi-dsi/aa/jaws/site.git
cp ./site/jaws_site/config.ini ../etc/jaws_site.ini
vi ../etc/jaws_site.ini
```

Fill in correct values provided by admin
```
chmod 400 ../etc/jaws_site.ini
cd site
pip install -r requirements.txt
```

---------------------------------

## Installing JAWS Client 
JAWS Client (Command-Line Interface)
NOTE: you only need to install JAWS Client on your laptop/desktop if you wish to process files accessible from your computer.  This configuration is intended for giving collaborators the ability to process their own data without granting them login access to servers.  Staff will generally use the JAWS Client installed on a server they SSH to.

Alternate installation methods exist but Anaconda is the easiest.
To install on your laptop:
Install (Globus Connect Personal)[https://www.globus.org/]
ensure it’s started after each reboot (on Mac: SysPrefs/User/Login Items)

Install (Anaconda)[https://www.anaconda.com/downloads]
python3 recommended
```
conda create -p .conda/envs/jaws python=3
cd ~/.conda/envs/jaws
mkdir -p opt etc/conda/activate.d
echo export JAWS_CLIENT_CONFIG=~/.conda/envs/jaws/etc/jaws_client.ini > etc/conda/activate.d/env_vars.sh
```
Install Cromwell’s womtool
```
conda activate jaws
conda install -c cyclus java-jdk
cd ~/.conda/envs/jaws/opt
mkdir cromwell
cd cromwell
```

Visit https://github.com/broadinstitute/cromwell/releases for latest release
```
wget https://github.com/broadinstitute/cromwell/releases/download/48/womtool-48.jar
wget https://github.com/broadinstitute/cromwell/releases/download/48/cromwell-48.jar
ln -s womtool-48.jar womtool.jar
ln -s cromwell-48.jar cromwell.jar
```

Install JAWS Client
```
cd ~/.conda/envs/jaws/opt
```

