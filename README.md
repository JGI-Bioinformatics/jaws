# cumentation for JGI Analysis Workflow Service (JAWS)

### Resources for JAWS Users
Full [documentation](https://jaws-docs.readthedocs.io) for running and installing JAWS is located here.

### Resources for JAWS Developers

### Installing JAWS site
See docs/install_Cli_and_Site.md  
See docs/startingServices.md  

### Installing all services 
You can install using two modes, developer mode and build mode:
you need to know when to you develop vs build.  
	develop => instant changes since it's just a symlink to your src
	build => the installed code will not change if src changes. This mode builds everything so that you will able to install the source using pip into your venv

if you're editing the code (i.e. developing) and wish to test your changes, it's inconvenient to reinstall via pip for every change, so you use develop.  when done and you're satisfied with the changes, use build to install a copy of your package.

For an install on your mac, everything is installed via develop.  So if you switch branches, all the servers will be running whatever version of the software is in that branch.  for shared installations, you'd want to use develop only while you're tweaking things, install via build before you leave because someone else may change the branch in the src dir.

Example commands using develop:

```
git clone https://code.jgi.doe.gov/advanced-analysis/jaws.git
module load python
python -m venv ~/venv/jaws-test  # do this once
source ~/venv/jaws-test/bin/activate
cd jaws/client && python setup.py develop
cd ../site && python setup.py develop
cd ../central && python setup.py develop
cd ../jtm && python setup.py develop
#To deactivate the venv:
deactivate
```

Example commands using build:

```
git clone https://code.jgi.doe.gov/advanced-analysis/jaws.git
module load python
python -m venv ~/venv/jaws-test  # do this once
source ~/venv/jaws-test/bin/activate
cd jaws/client && python setup.py build && pip install .
cd ../site && python setup.py build && pip install .
cd ../central && python setup.py build && pip install .
cd ../jtm && python setup.py build && pip install .

#To deactivate the venv:
deactivate
```

### cromwell-utils
Example installing cromwell-utils using build mode

```
# do this once
module load python/3.7
python -m venv ~/cromvenv
cd <path_to_repo>/jaws/cromwell_utilities/
python setup.py [develop|build]
pip install .  # run this if you used build in last step

# now do this every time
export TMPDIR=/"global/scratch/$USER"
export CROMWELL_URL=localhost:50010
source ~/cromvenv/bin/activate
cromwell-utils
```

## Contributing
* Nicole Free <nlfree@lbl.gov>
* Jeff Froula <jlfroula@lbl.gov>  
* Edward Kirton <eskirton@lbl.gov>  
* Mario Melara <mamelara@lbl.gov>  
* Georg Rath <gbrath@lbl.gov>  
* Hugh Salamon <hsalamon@lbl.gov> 
* Seung-Jin Sul <ssul@lbl.gov>   
* Stephan Trong <strong@lbl.gov>  
