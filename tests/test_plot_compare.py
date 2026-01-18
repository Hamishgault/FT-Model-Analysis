import matplotlib
matplotlib.use("Agg")
import numpy as np
from pathlib import Path

from kinetics import BrubachModel
from utils.io import load_experimental_json
from utils.parameters import REACTOR
from utils.plotting import plot_model_vs_experiment
from reactor.pfr import PFR
from utils.compare import make_inlet_from_experiment


def test_plot_overlay_generates_file(tmp_path: Path):
    records = load_experimental_json("data/raw/sample_experiments.json")
    rec = records[0]
    model = BrubachModel()

    Fi = make_inlet_from_experiment(rec)
    pfr = PFR(REACTOR["A"], REACTOR["L"], model, nu=None)
    z_eval = np.linspace(0.0, REACTOR["L"], 101)
    sol = pfr.run(rec.get("T", 523.15), 1e5, Fi, z_eval=z_eval)
    z, y = sol.t, sol.y

    out = tmp_path / "overlay.png"
    fig = plot_model_vs_experiment(z, y, rec, {
        "CO": 0, "H2": 1, "CH4": 2, "C2_4": 3, "C5plus": 4, "H2O": 5, "CO2": 6
    }, str(out))

    assert out.exists()
    assert out.stat().st_size > 0
