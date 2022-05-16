## Repository for JAWS docs
This content was created by using sphinx and hosted by readthedocs

### Install sphinx

```
apt-get install python3-sphinx ## Linux
brew install sphinx-doc ## macOS
```

### I created conda env

```
conda create --name sphinx
conda activate sphinx
pip install sphinx
pip install recommonmark
pip install -U sphinxcontrib-confluencebuilder
pip install sphinx_rtd_theme 
```
  
### creating the docs
follow sphinx tutorial in readthedocs.org
https://docs.readthedocs.io/en/latest/intro/getting-started-with-sphinx.html

Sphinx uses reStructuredText:
http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

### create initial directory structure
```
sphinx-quickstart
```

Now start adding content to the pages 

### build html from source (rst or md) files
```
make html
```

### and quick view
```
open build/html/index.html
```

Once you are satisfied, add changes to this gitlab repo


## Register gitlab with readthedocs
https://docs.readthedocs.io/en/latest/intro/getting-started-with-sphinx.html

## view published docs
https://jaws-docs.readthedocs.io/en/latest/

# Publish pages to confluence

### followed tutorial
https://sphinxcontrib-confluencebuilder.readthedocs.io/en/latest/tutorial.html


### install confluence builder
```
on cori
# working dir
/global/cscratch1/sd/jfroula/

# create env
/global/dna/shared/data/jfroula/miniconda3/bin/virtualenv testenv
. ./testenv/bin/activate
cd testenv/

# install dependencies
pip install sphinxcontrib-confluencebuilder
pip install recommonmark

# clone repo

# add confluence stuff to conf.py

# this command should add pages to confluence (space JAWS)
make confluence
```
