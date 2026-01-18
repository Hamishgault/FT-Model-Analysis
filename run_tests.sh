#!/usr/bin/env bash
set -e

# activate local .venv and run tests
. .venv/bin/activate
pytest --cov=./ -q
