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

## Continuous integration (GitHub Actions) 🤖

A GitHub Actions workflow is provided in `.github/workflows/python-tests.yml` that runs the test suite on push and pull requests for Python 3.11 and 3.12.

Per-model stoichiometry

- Models may expose a `nu` attribute (numpy array shape `n_species x n_rxns`). When present, the `PFR` will use `model.nu` in preference to the global `nu` passed to the reactor. This enables each kinetic module to declare its own stoichiometric mapping cleanly.

Experimental data and comparisons 📁

- Add raw experimental files to `data/raw/` and a per-dataset metadata JSON to `data/metadata/`.
- Use `utils.io.load_experimental_json` and `utils.io.process_sample_records` to load and preprocess datasets into a standard table.
- Use `utils.compare.compare_dataset` to run a model against the dataset and compute simple comparison metrics (CO outlet, selectivity differences).

You can run tests locally using the provided helper scripts:

PowerShell (Windows):

```powershell
.\run_tests.ps1
```

Bash (macOS / Linux):

```bash
./run_tests.sh
```

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
