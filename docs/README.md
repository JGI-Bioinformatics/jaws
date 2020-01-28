## Repository for JAWS docs
This content was created by using sphinx and hosted by readthedocs

### I created conda env

```
conda create sphinx
conda activate sphinx
pip install -y sphinx
pip install recommonmark
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


### Register gitlab with readthedocs
https://docs.readthedocs.io/en/latest/intro/getting-started-with-sphinx.html

### view published docs
https://jaws-sphinx.readthedocs.io/en/latest/

### Add new repo to official JAWS repository
`jgi-workflows/docs`
