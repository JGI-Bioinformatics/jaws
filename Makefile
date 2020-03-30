VERSION := $(shell git describe --always --tags --abbrev=6 --dirty="-dev")
Q := $(if $V,,@)

## Package Section BEGIN
.PHONY: pkg-requirements
pkg-requirements:
	$(if $(shell which wheel),,$(error "Packaging needs Python wheel installed. Please run 'pip install wheel'"))

.PHONY: pkg-site
pkg-site: pkg-requirements
	$Q cd site && python setup.py bdist_wheel

.PHONY: pkg-central
pkg-central: pkg-requirements
	$Q cd central && python setup.py bdist_wheel

.PHONY: pkg-client
pkg-client: pkg-requirements
	$Q cd client && python setup.py bdist_wheel

.PHONY: pkg-jtm
pkg-jtm: pkg-requirements
	$Q cd jtm && python setup.py bdist_wheel

.PHONY: pkg
pkg: pkg-site pkg-central pkg-client pkg-jtm
## Package Section END

## Test Section BEGIN
.PHONY: test-requirements
test-requirements:
	$(if $(shell which flake8),,$(error "Testing needs flake8 installed. Please run 'pip install flake8'"))
	$(if $(shell which pytest),,$(error "Testing needs pytest installed. Please run 'pip install pytest'"))

.PHONY: test-site
test-site: test-requirements
	$Q flake8 site

.PHONY: test-central
test-central: test-requirements
	$Q flake8 central

.PHONY: test-client
test-client: test-requirements
	$Q flake8 client
	$Q cd client && pytest

.PHONY: test-jtm
test-jtm: test-requirements
	$Q flake8 jtm

.PHONY: test
test: test-site test-central test-client test-jtm
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
	$Q echo "open docs/build/html/index.html"
## Doc Section END
