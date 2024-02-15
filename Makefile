VERSION := $(shell git describe --always --tags --abbrev=0)
Q := $(if $V,,@)

init:
	pip install --no-cache-dir -r requirements.txt
	pip install --no-cache-dir -e .

init-dev: init
	pip install .[dev]

update-deps:
	pip install --upgrade pip-tools pip setuptools
	pip-compile --upgrade --build-isolation \
		--allow-unsafe --resolver=backtracking --strip-extras \
		--output-file requirements.txt \
		pyproject.toml

update: update-deps init


up-dev:
	podman-compose up --build --force-recreate --detach --remove-orphans


down-dev:
	podman-compose down
## Package Section BEGIN

.PHONY: pkg-requirements
pkg-requirements:
	$(if $(shell which wheel),,$(error "Packaging needs Python wheel installed. Please run 'pip install wheel'"))

.PHONY: pkg-rpc
pkg-rpc: pkg-requirements
        $Q pip install -r requirements.txt
	$Q python -m build

.PHONY: pkg-site
pkg-site: pkg-requirements
        $Q pip install -r requirements.txt
	$Q python -m build

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
	$Q flake8 src/jaws_rpc
	$Q python -m pytest --cov=jaws_rpc --junitxml=rpc.xml rpc/tests/ && coverage xml

.PHONY: test-site
test-site: test-requirements
	$Q flake8 src/jaws_site
	$Q python -m pytest --cov=jaws_site --junitxml=site.xml site/tests/ && coverage xml

.PHONY: test
test: test-rpc test-site
## Test Section END
