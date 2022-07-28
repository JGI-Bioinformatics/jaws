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

.PHONY: pkg-poetry
pkg-poetry: pkg-rpc-poetry pkg-site-poetry

.PHONY: pkg-rpc
pkg-rpc: pkg-requirements
	$Q cd rpc && python setup.py bdist_wheel

.PHONY: pkg-site
pkg-site: pkg-requirements
	$Q cd rpc && python setup.py bdist_wheel && cd ../site && python setup.py bdist_wheel

.PHONY: pkg
pkg: pkg-rpc pkg-site

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

.PHONY: test-condor
test-condor: test-requirements
	$Q flake8 condor
	$Q cd condor && python -m pytest --cov=jaws_condor --junitxml=condor.xml tests/ && coverage xml

.PHONY: test
test: test-rpc test-site test-condor
## Test Section END
