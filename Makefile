VERSION := $(shell git describe --always --tags --abbrev=0)
Q := $(if $V,,@)

init:
	pip install -r site/requirements.txt


update-deps:
	pip install --upgrade pip-tools pip setuptools
	pip-compile --upgrade --build-isolation \
		--allow-unsafe --resolver=backtracking --strip-extras \
		--output-file site/requirements.txt \
		site/pyproject.toml rpc/pyproject.toml pubsub/pyproject.toml


update: update-deps init


up-dev:
	docker compose up --build --force-recreate --detach --remove-orphans


down-dev:
	docker compose down


.PHONY: update-deps init update
## Package Section BEGIN
.PHONY: pkg-requirements
pkg-requirements:
	$(if $(shell which wheel),,$(error "Packaging needs Python wheel installed. Please run 'pip install wheel'"))

.PHONY: pkg-rpc
pkg-rpc: pkg-requirements
	$Q cd rpc && python -m build

.PHONY: pkg-site
pkg-site: pkg-requirements
	$Q cd rpc && python -m build && cd ../site && python -m build

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
	$Q cd site && python -m pytest --cov=jaws_site --junitxml=site.xml tests/ && coverage xml

.PHONY: test
test: test-rpc test-site
## Test Section END
