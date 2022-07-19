VERSION := $(shell git describe --always --tags --abbrev=0)
Q := $(if $V,,@)

## Package Section BEGIN
.PHONY: pkg-requirements
pkg-requirements:
	$(if $(shell which wheel),,$(error "Packaging needs Python wheel installed. Please run 'pip install wheel'"))

.PHONY: pkg-poetry-requirements
pkg-poetry-requirements:
	$(if $(shell which poetry),,$(error "Packaging needs Python poetry module installed. Please run 'pip install poetry'"))

.PHONY: pkg-rpc-poetry
pkg-rpc-poetry: pkg-poetry-requirements
	$Q cd rpc && poetry version $(VERSION) && poetry build

.PHONY: pkg-site-poetry
pkg-site-poetry: pkg-poetry-requirements
	$Q cd site && poetry version $(VERSION) && poetry build 

.PHONY: pkg-central-poetry
pkg-central-poetry: pkg-poetry-requirements
	$Q cd central && poetry version $(VERSION) && poetry build

.PHONY: pkg-jtm-poetry
pkg-jtm-poetry: pkg-poetry-requirements
	$Q cd jtm && poetry version $(VERSION) && poetry build

.PHONY: pkg-poetry
pkg-poetry: pkg-rpc-poetry pkg-site-poetry pkg-central-poetry pkg-jtm-poetry

.PHONY: pkg-rpc
pkg-rpc: pkg-requirements
	$Q cd rpc && python setup.py bdist_wheel

.PHONY: pkg-site
pkg-site: pkg-requirements
	$Q cd rpc && python setup.py bdist_wheel && cd ../site && python setup.py bdist_wheel

.PHONY: pkg-central
pkg-central: pkg-requirements
	$Q cd rpc && python setup.py bdist_wheel && cd ../central && python setup.py bdist_wheel

.PHONY: pkg-jtm
pkg-jtm: pkg-requirements
	$Q cd rpc && python setup.py bdist_wheel && cd ../jtm && python setup.py bdist_wheel

.PHONY: pkg-parsl
pkg-parsl: pkg-requirements
	$Q cd rpc && python setup.py bdist_wheel && cd ../parsl && python setup.py bdist_wheel

.PHONY: pkg
pkg: pkg-rpc pkg-site pkg-central pkg-jtm pkg-parsl

## Package Section END

## Test Section BEGIN
.PHONY: test-requirements
test-requirements:
	$(if $(shell which flake8),,$(error "Testing needs flake8 installed. Please run 'pip install flake8'"))
	$(if $(shell which pytest),,$(error "Testing needs pytest installed. Please run 'pip install pytest'"))

.PHONY: test-rpc
test-rpc: test-requirements
	$Q flake8 rpc
	$Q cd rpc && python -m pytest --cov=jaws_rpc --junitxml=rpc.xml tests/ && coverage xml

.PHONY: test-site
test-site: test-requirements
	$Q flake8 site
	$Q cd site && python -m pytest --cov=jaws_site --cov=jaws_site/datatransfer_plugins --junitxml=site.xml tests/ && coverage xml

.PHONY: test-central
test-central: test-requirements
	$Q flake8 central
	$Q cd central && python -m pytest --cov=jaws_central --junitxml=central.xml tests/ && coverage xml

.PHONY: test-jtm
test-jtm: test-requirements
	$Q flake8 jtm
	$Q cd jtm && python -m pytest --cov=jaws_jtm --junitxml=jtm.xml tests/ && coverage xml

.PHONY: test-condor
test-condor: test-requirements
	$Q flake8 condor
	$Q cd condor && python -m pytest --cov=jaws_condor --junitxml=condor.xml tests/ && coverage xml

.PHONY: test-parsl
test-parsl: test-requirements
	$Q flake8 parsl
	$Q cd parsl && python -m pytest --cov=jaws_parsl --junitxml=parsl.xml tests/ && coverage xml

.PHONY: test
test: test-rpc test-site test-central test-jtm test-parsl test-condor
## Test Section END

## Doc Section BEGIN
.PHONY: doc-requirements
doc-requirements:
	$(if $(shell python -c 'import sphinx; print("ok")' 2>/dev/null),,$(error "Testing needs sphinx installed. Please run 'pip install sphinx'"))
	$(if $(shell python -c 'from sphinxcontrib import confluencebuilder; print("ok")' 2>/dev/null),,$(error "Building sphinx docs needs sphinxcontrib-confluencebuilder installed. Please run 'pip install sphinxcontrib-confluencebuilder'"))
	$(if $(shell python -c 'import recommonmark; print("ok")' 2>/dev/null),,$(error "Building sphinx docs needs recommonmark installed. Please run 'pip install recommonmark'"))
	$(if $(shell python -c 'import sphinx_rtd_theme; print("ok")' 2>/dev/null),,$(error "Building sphinx docs needs sphinx_rtd_theme installed. Please run 'pip install sphinx_rtd_theme'"))

.PHONY: docs
docs: doc-requirements
	$Q cd docs/sphinx && make html
	$Q echo "open docs/sphinx/build/html/index.html"
## Doc Section END
