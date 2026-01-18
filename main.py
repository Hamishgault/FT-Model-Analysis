"""CLI-enabled runner for FT model PFR simulations.

This script provides a small command-line interface to run one or more kinetic
models in the isothermal PFR and inspect/save simple outputs (conversion
profiles and selectivity summary).

Usage (examples):
    python main.py                       # runs all models and shows plots
    python main.py --models PowerLaw     # run only PowerLaw model
    python main.py --save --outdir out   # save plots to `out/`
"""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import matplotlib.pyplot as plt

from kinetics import (
    PowerLawModel,
    IglesiaCOInsertionModel,
    SteynbergCarbideModel,
    VanDerLaanAlkenylModel,
)
from reactor.pfr import PFR
from utils.parameters import REACTOR, F_INLET, T0, P0, MODEL_PARAMS
from utils.plotting import plot_conversion, plot_selectivity_bar
from utils.species import SPECIES_IDX, NU, compute_selectivity


AVAILABLE_MODELS = {
    "PowerLaw": PowerLawModel,
    "Iglesia": IglesiaCOInsertionModel,
    "Steynberg": SteynbergCarbideModel,
    "VanDerLaan": VanDerLaanAlkenylModel,
}


def run_model(model, nu, reactor_params, F0, T, P, z_points=201):
    pfr = PFR(reactor_params["A"], reactor_params["L"], model, nu)
    z_eval = np.linspace(0.0, reactor_params["L"], z_points)
    sol = pfr.run(T, P, F0, z_eval=z_eval)
    return sol.t, sol.y


def summary_from_profiles(F0, y_out):
    F_out = y_out[:, -1]
    CO_conv = (F0[SPECIES_IDX["CO"]] - F_out[SPECIES_IDX["CO"]]) / max(F0[SPECIES_IDX["CO"]], 1e-12)
    sel = compute_selectivity(F0, F_out)
    return {"CO_outlet": float(F_out[SPECIES_IDX["CO"]]), "CO_conversion": float(CO_conv), "selectivity": sel}


def run_and_report(model_names: Iterable[str], save: bool = False, outdir: Path | None = None, z_points: int = 201):
    models = {name: AVAILABLE_MODELS[name](MODEL_PARAMS.get(name.lower(), None)) for name in model_names}

    results = {}
    summaries = {}

    for name, model in models.items():
        logging.info("Running model: %s", name)
        z, y = run_model(model, NU, REACTOR, F_INLET.copy(), T0, P0, z_points=z_points)
        results[name] = (z, y)
        summaries[name] = summary_from_profiles(F_INLET, y)
        logging.info("Model %s summary: %s", name, summaries[name])

    # Plot conversion (CO)
    plt.figure()
    for name, (z, y) in results.items():
        F0_co = F_INLET[SPECIES_IDX["CO"]]
        conv = (F0_co - y[SPECIES_IDX["CO"], :]) / max(F0_co, 1e-12)
        plt.plot(z, conv, label=name)
    plt.xlabel("Reactor length, m")
    plt.ylabel("CO conversion")
    plt.legend()
    plt.grid(True)
    plt.title("CO conversion vs reactor length")

    if save and outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        conv_path = outdir / "conversion_co.png"
        plt.savefig(conv_path, dpi=200)
        logging.info("Saved conversion plot to %s", conv_path)

    # Selectivity bar chart
    plt.figure()
    plot_selectivity_bar({name: summaries[name]["selectivity"] for name in summaries})
    if save and outdir is not None:
        sel_path = outdir / "selectivity.png"
        plt.savefig(sel_path, dpi=200)
        logging.info("Saved selectivity plot to %s", sel_path)

    # Save JSON summary
    if save and outdir is not None:
        summary_path = outdir / "summaries.json"
        with open(summary_path, "w") as fh:
            json.dump(summaries, fh, indent=2)
        logging.info("Saved summaries to %s", summary_path)

    # Show plots
    plt.show()
    return results, summaries


def parse_args():
    parser = argparse.ArgumentParser(description="Run FT kinetic models in a 1D PFR and compare outputs")
    parser.add_argument("--models", type=str, default="all", help="Comma-separated model names to run (default: all)")
    parser.add_argument("--save", action="store_true", help="Save plots and summaries to disk")
    parser.add_argument("--outdir", type=str, default="out", help="Output directory when saving plots/summaries")
    parser.add_argument("--zpoints", type=int, default=201, help="Number of axial points to evaluate")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s: %(message)s")

    if args.models.lower() == "all":
        model_names = list(AVAILABLE_MODELS.keys())
    else:
        model_names = [m.strip() for m in args.models.split(",") if m.strip() in AVAILABLE_MODELS]
        if not model_names:
            logging.error("No valid model names provided. Available: %s", ", ".join(AVAILABLE_MODELS.keys()))
            return

    outdir = Path(args.outdir) if args.save else None
    results, summaries = run_and_report(model_names, save=args.save, outdir=outdir, z_points=args.zpoints)

    # Print brief summary
    for name, s in summaries.items():
        print(f"Model: {name}")
        print(f"  CO outlet: {s['CO_outlet']:.4f} mol/s  CO conversion: {s['CO_conversion']*100:.2f}%")
        print(f"  Selectivity: {s['selectivity']}\n")


if __name__ == "__main__":
    main()