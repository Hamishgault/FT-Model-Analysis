# FT Model Analysis

A modular Python framework to compare Fischer–Tropsch (FT) kinetic models in a 1D plug-flow reactor (PFR).

## Quick start ✅

Recommended: use the provided helper scripts to create a `.venv` and install dependencies.

PowerShell (Windows):

```powershell
# create venv and install
.\setup_env.ps1
# activate
. .venv\Scripts\Activate.ps1
python main.py
```

Bash (macOS / Linux):

```bash
./setup_env.sh
source .venv/bin/activate
python main.py
```

Manual:

```bash
python -m venv .venv
.venv\Scripts\python.exe -m pip install --upgrade pip  # Windows
.venv\Scripts\pip.exe install -r requirements.txt      # Windows
source .venv/bin/activate                              # Unix
python main.py
```

## Project layout

- `kinetics/` - kinetic model interface and model modules
- `reactor/` - PFR solver
- `utils/` - species, parameters, plotting helpers
- `main.py` - example driver

## Notes

- Current kinetics are placeholders; replace rate laws with literature expressions as needed.
- For best reproducibility, run inside the `.venv` created by the helper scripts.

## Basic workflow (CLI) ▶️

The project includes a CLI-enabled `main.py` you can run to execute one or more models and visualize/save results.

Examples:

Run all models and show plots:

```bash
python main.py
```

Run a single model and save outputs to `out/`:

```bash
python main.py --models PowerLaw --save --outdir out
```

Run a comma-separated list of models:

```bash
python main.py --models PowerLaw,Iglesia --save
```

By default, plots are shown interactively. Use `--save` to also write PNGs and a `summaries.json` into `--outdir`.

## Tests ✅

Run unit tests using `pytest` from the project root (inside the activated `.venv`):

PowerShell (Windows):

```powershell
. .venv\Scripts\Activate.ps1
pytest -q
```

Bash (macOS/Linux):

```bash
source .venv/bin/activate
pytest -q
```

To run with coverage and produce a simple report:

```bash
pytest --cov=./ -q
```
