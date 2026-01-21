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
    BrubachModel,
)
from reactor.pfr import PFR
from utils.parameters import REACTOR, F_INLET, T0, P0, MODEL_PARAMS
from utils.plotting import plot_conversion, plot_selectivity_bar
from utils.species import SPECIES_IDX, NU, compute_selectivity


AVAILABLE_MODELS = {
    "Brubach": BrubachModel,
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


def run_and_report(model_names: Iterable[str], save: bool = False, outdir: Path | None = None, z_points: int = 201, diagnostics: bool = False, diag_outdir: str | None = None):
    models = {name: AVAILABLE_MODELS[name](MODEL_PARAMS.get(name.lower(), None)) for name in model_names}

    results = {}
    summaries = {}

    for name, model in models.items():
        # Apply temporary parameter overrides for Brubach if provided
        if name.lower() == "brubach":
            overrides = MODEL_PARAMS.get("brubach", {}) or {}
            # Respect any CLI overrides added to overrides dict earlier
            overrides.update({k: v for k, v in MODEL_PARAMS.get("brubach_overrides", {}).items() if v is not None})
            if overrides:
                model = AVAILABLE_MODELS[name](overrides)

        logging.info("Running model: %s", name)
        z, y = run_model(model, NU, REACTOR, F_INLET.copy(), T0, P0, z_points=z_points)
        results[name] = (z, y)
        summaries[name] = summary_from_profiles(F_INLET, y)
        logging.info("Model %s summary: %s", name, summaries[name])

        # if diagnostics requested, compute per-z diagnostics and optionally save plots
        if diagnostics:
            try:
                from utils.plotting import compute_model_diagnostics, plot_diagnostics
                diag = compute_model_diagnostics(model, T0, P0, z, y)
                print(f"Diagnostics for model {name} at z points [0, L/2, L]:")
                for idx, zpt in enumerate((0, len(z) // 2, -1)):
                    print(f" z={z[zpt]:.3f} m: theta_CH2={diag['theta_CH2'][zpt]:.3e}, r_growth={diag['r_growth'][zpt]:.3e}, r_ch4={diag['r_ch4'][zpt]:.3e}, r_c2_4={diag['r_c2_4'][zpt]:.3e}, r_c5p={diag['r_c5p'][zpt]:.3e}")

                dout = None
                if diag_outdir is not None:
                    dout = Path(diag_outdir)
                    dout.mkdir(parents=True, exist_ok=True)
                outpath = dout / f"{name}_diagnostics.png" if dout is not None else None
                plot_diagnostics(z, diag, outpath)
                logging.info("Saved diagnostics plot to %s", outpath) if outpath is not None else None
            except Exception as e:
                logging.warning("Failed to compute diagnostics for %s: %s", name, e)

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
    parser.add_argument("--compare", type=str, default=None, help="Path to experimental dataset (JSON or CSV) to compare against")
    parser.add_argument("--menu", action="store_true", help="Start interactive menu to access common functions")
    # Temporary parameter overrides for Brubach tuning
    parser.add_argument("--k6", type=float, default=None, help="Override k6 (CO hydrogenation RDS) in BrubachParams")
    parser.add_argument("--k8", type=float, default=None, help="Override k8 (chain growth) in BrubachParams")
    parser.add_argument("--k9a", type=float, default=None, help="Override k9a (methane termination) in BrubachParams")
    parser.add_argument("--gamma10", type=float, default=None, help="Override Gamma10 (J/mol) in BrubachParams")
    parser.add_argument("--cat-loading", type=float, default=None, help="Override catalyst loading (g/m^3) in BrubachParams")
    parser.add_argument("--diagnostics", action="store_true", help="Compute and plot diagnostic variables (theta_CH2, rates) vs reactor length")
    parser.add_argument("--diag-outdir", type=str, default=None, help="Directory to save diagnostic plots (if --diagnostics)")
    return parser.parse_args()


def interactive_menu():
    """Simple interactive menu to access common functions.

    Options:
    1) Run models
    2) Compare with dataset (JSON/CSV)
    3) Run tests (pytest)
    4) Exit
    """
    import subprocess

    while True:
        print("\n=== FT Model Analysis - Interactive Menu ===")
        print("1) Run models")
        print("2) Compare with experimental dataset (JSON/CSV)")
        print("3) Run tests (pytest)")
        print("4) Exit")
        try:
            choice = input("Select option [1-4]: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nExiting menu.")
            return

        if choice == "1":
            models = input("Enter models to run (comma-separated or 'all') [all]: ").strip() or "all"
            save_ans = input("Save outputs? [y/N]: ").strip().lower()
            outdir = input("Output directory [out]: ").strip() or "out"
            zpoints = input("Axial points (int) [201]: ").strip() or "201"
            try:
                zpoints = int(zpoints)
            except ValueError:
                print("Invalid zpoints, using 201")
                zpoints = 201
            if models.lower() == "all":
                model_names = list(AVAILABLE_MODELS.keys())
            else:
                model_names = [m.strip() for m in models.split(",") if m.strip() in AVAILABLE_MODELS]
                if not model_names:
                    print(f"No valid models provided. Available: {', '.join(AVAILABLE_MODELS.keys())}")
                    continue
            run_and_report(model_names, save=(save_ans == "y"), outdir=Path(outdir), z_points=zpoints)

        elif choice == "2":
            path = input("Path to experimental dataset (JSON or CSV): ").strip()
            if not path:
                print("No path provided. Aborting.")
                continue
            comp_outdir = input("Output dir for comparisons [out/comparisons]: ").strip() or "out/comparisons"
            try:
                from utils.io import load_experimental_json, load_experimental_csv
                from utils.compare import make_inlet_from_experiment, compare_dataset, save_comparison_results
                from utils.plotting import plot_model_vs_experiment
            except Exception as e:
                print(f"Failed to import comparison utilities: {e}")
                continue

            cmpath = Path(path)
            if cmpath.suffix.lower() == ".json":
                try:
                    records = load_experimental_json(str(cmpath))
                except Exception as e:
                    print(f"Failed to load JSON: {e}")
                    continue
            else:
                try:
                    records = load_experimental_csv(str(cmpath))
                except Exception as e:
                    print(f"Failed to load CSV: {e}")
                    continue

            comp_outdir = Path(comp_outdir)
            comp_outdir.mkdir(parents=True, exist_ok=True)

            # Run models for each record and save overlays
            for name, model in {n: AVAILABLE_MODELS[n](MODEL_PARAMS.get(n.lower(), None)) for n in AVAILABLE_MODELS}.items():
                for rec in records:
                    Fi = make_inlet_from_experiment(rec)
                    pfr = PFR(REACTOR["A"], REACTOR["L"], model, nu=None)
                    sol = pfr.run(rec.get("T", T0), P0, Fi, z_eval=np.linspace(0.0, REACTOR["L"], 201))
                    outpath = comp_outdir / f"{name}_{rec.get('id')}_overlay.png"
                    plot_model_vs_experiment(sol.t, sol.y, rec, SPECIES_IDX, str(outpath))
            # aggregate & save comparison table
            try:
                comp_df = compare_dataset(records, AVAILABLE_MODELS[list(AVAILABLE_MODELS.keys())[0]](MODEL_PARAMS.get(list(AVAILABLE_MODELS.keys())[0].lower(), None)), REACTOR)
                save_comparison_results(comp_df, str(comp_outdir / "comparison_summary.csv"))
                print(f"Saved comparison overlays and summary to {comp_outdir}")
            except Exception as e:
                print(f"Failed to aggregate comparison table: {e}")

        elif choice == "3":
            print("Running pytest...")
            try:
                subprocess.run(["pytest", "-q"], check=False)
            except Exception as e:
                print(f"Failed to run pytest: {e}")

        elif choice == "4":
            print("Exiting menu.")
            return
        else:
            print("Invalid option, please try again.")


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

    # Apply CLI overrides for Brubach parameters if provided
    if args.k6 is not None or args.k8 is not None or args.k9a is not None or args.gamma10 is not None or args.cat_loading is not None:
        MODEL_PARAMS.setdefault("brubach_overrides", {})
        if args.k6 is not None:
            MODEL_PARAMS["brubach_overrides"]["k6"] = args.k6
        if args.k8 is not None:
            MODEL_PARAMS["brubach_overrides"]["k8"] = args.k8
        if args.k9a is not None:
            MODEL_PARAMS["brubach_overrides"]["k9a"] = args.k9a
        if args.gamma10 is not None:
            MODEL_PARAMS["brubach_overrides"]["Gamma10"] = args.gamma10
        if args.cat_loading is not None:
            MODEL_PARAMS["brubach_overrides"]["cat_loading"] = args.cat_loading

    # If interactive menu requested, start it
    if args.menu:
        interactive_menu()
        return

    diagnostics = args.diagnostics
    diag_outdir = args.diag_outdir

    results, summaries = run_and_report(model_names, save=args.save, outdir=outdir, z_points=args.zpoints, diagnostics=diagnostics, diag_outdir=diag_outdir)

    # If user asked to compare with experimental dataset, run comparisons and save overlay plots
    if args.compare:
        from utils.io import load_experimental_json, load_experimental_csv
        from utils.plotting import plot_model_vs_experiment
        from utils.compare import compare_dataset, save_comparison_results

        comp_outdir = Path(args.outdir or "out") / "comparisons"
        comp_outdir.mkdir(parents=True, exist_ok=True)

        # load records from JSON or CSV
        cmpath = Path(args.compare)
        if cmpath.suffix.lower() == ".json":
            records = load_experimental_json(str(cmpath))
        else:
            records = load_experimental_csv(str(cmpath))

        # run comparisons for each model and plot overlays
        for name, model in {n: AVAILABLE_MODELS[n](MODEL_PARAMS.get(n.lower(), None)) for n in model_names}.items():
            for rec in records:
                # build inlet and run short PFR for plotting
                from utils.compare import make_inlet_from_experiment
                Fi = make_inlet_from_experiment(rec)
                pfr = PFR(REACTOR["A"], REACTOR["L"], model, nu=None)
                z, y = pfr.run(rec.get("T", T0), P0, Fi, z_eval=np.linspace(0.0, REACTOR["L"], args.zpoints)).t, pfr.run(rec.get("T", T0), P0, Fi, z_eval=np.linspace(0.0, REACTOR["L"], args.zpoints)).y
                # save overlay
                outpath = comp_outdir / f"{name}_{rec.get('id')}_overlay.png"
                plot_model_vs_experiment(z, y, rec, SPECIES_IDX, str(outpath))

        # aggregate & save comparison table
        comp_df = compare_dataset(records, AVAILABLE_MODELS[model_names[0]](MODEL_PARAMS.get(model_names[0].lower(), None)), REACTOR, z_points=args.zpoints)
        save_comparison_results(comp_df, str(comp_outdir / "comparison_summary.csv"))
        print(f"Saved comparison overlays and summary to {comp_outdir}")

    # Print brief summary
    for name, s in summaries.items():
        print(f"Model: {name}")
        print(f"  CO outlet: {s['CO_outlet']:.4f} mol/s  CO conversion: {s['CO_conversion']*100:.2f}%")
        print(f"  Selectivity: {s['selectivity']}\n")


if __name__ == "__main__":
    main()