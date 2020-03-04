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

.PHONY: pkg-auth
pkg-auth: pkg-requirements
	$Q cd auth && python setup.py bdist_wheel

.PHONY: pkg
pkg: pkg-site pkg-central pkg-client
## Package Section END

## Test Section BEGIN
.PHONY: test-requirements
test-requirements:
	$(if $(shell which flake8),,$(error "Testing needs flake8 installed. Please run 'pip install flake8'"))

.PHONY: test-site
test-site: test-requirements
	$Q flake8 site

.PHONY: test-central
test-central: test-requirements
	$Q flake8 central

.PHONY: test-client
test-client: test-requirements
	$Q flake8 client

.PHONY: test-auth
test-auth: test-requirements
	$Q flake8 auth
	
.PHONY: test
test: test-site test-central test-client test-auth
## Test Section END
