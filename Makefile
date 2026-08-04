.PHONY: install install-all test lint typecheck validate publication-check figures build ci

install:
	python -m pip install -e ".[dev,xenium]"

install-all:
	python -m pip install -e ".[dev,xenium,viz]"

test:
	pytest -q

lint:
	ruff check src tests scripts
	ruff format --check src tests scripts

typecheck:
	mypy src/reach_gap

validate:
	reach validate-v080-results results/evidence_synthesis_v0.8

publication-check:
	python scripts/validate_github_publication.py

figures:
	python scripts/generate_portfolio_figures.py

build:
	python -m build

ci: lint typecheck test validate publication-check
