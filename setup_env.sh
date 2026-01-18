#!/usr/bin/env bash
set -e
PYTHON=${1:-python3}

echo "Creating virtual environment at .venv..."
$PYTHON -m venv .venv

echo "Activating .venv and installing requirements..."
. .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt

echo "Done. Activate with: source .venv/bin/activate"
