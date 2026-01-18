"""Plotting helpers for FT model comparisons.

These are simple wrappers around matplotlib to produce conversion and
selectivity plots. Kept small so they can be adapted or replaced.
"""
from typing import Dict, Iterable
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def plot_conversion(z: np.ndarray, F_profiles: Dict[str, np.ndarray], species_idx: int, label_map: dict = None):
    """Plot conversion of a species vs reactor length for multiple models.

    Parameters
    ----------
    z : np.ndarray
        Reactor axial positions
    F_profiles : dict
        Mapping model name -> array of species flows (n_species, n_z)
    species_idx : int
        Index of species to compute conversion for
    """
    plt.figure()
    for name, F in F_profiles.items():
        F = np.atleast_2d(F)
        F0 = F[species_idx, 0]
        conv = (F0 - F[species_idx, :]) / max(F0, 1e-12)
        plt.plot(z, conv, label=name)
    plt.xlabel("Reactor length, m")
    plt.ylabel("Conversion")
    plt.legend()
    plt.grid(True)
    plt.title("Conversion vs Reactor Length")


def plot_selectivity_bar(selectivities: Dict[str, Dict[str, float]]):
    """Bar chart of selectivity for each model (final outlet)."""
    labels = list(selectivities.keys())
    groups = ["CH4", "C2_4", "C5+"]

    x = np.arange(len(labels))
    width = 0.2

    fig, ax = plt.subplots()
    for i, g in enumerate(groups):
        vals = [selectivities[name].get(g, 0.0) for name in labels]
        ax.bar(x + (i - 1) * width, vals, width, label=g)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel("Selectivity")
    ax.set_title("Product Selectivity (outlet)")
    ax.legend()
    plt.tight_layout()


def plot_model_vs_experiment(z_model, F_model, exp_record: dict, species_keys: dict, outpath: str | None = None):
    """Overlay model profiles and experimental points for a single record.

    Parameters
    - z_model: axial positions (1D array)
    - F_model: array of shape (n_species, n_z)
    - exp_record: processed experimental record/dict with keys like 'CO_out', 'sel_CH4' etc.
    - species_keys: mapping species name -> index in F_model (e.g., SPECIES_IDX)
    - outpath: optional path to save figure (PNG). If None, figure is shown.
    """
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    # 1) CO conversion profile (model line, experimental point as marker)
    iCO = species_keys["CO"]
    F0_co = F_model[iCO, 0]
    conv = (F0_co - F_model[iCO, :]) / max(F0_co, 1e-12)
    axes[0].plot(z_model, conv, label="Model")
    # experimental CO_out is outlet fraction; plot as a single point at max(z)
    if exp_record.get("CO_out") is not None:
        co_out = exp_record["CO_out"]
        conv_exp = (F0_co - co_out) / max(F0_co, 1e-12)
        axes[0].plot([z_model[-1]], [conv_exp], marker="o", color="C1", label="Exp")
    axes[0].set_xlabel("Reactor length, m")
    axes[0].set_ylabel("CO conversion")
    axes[0].legend()
    axes[0].grid(True)

    # 2) Selectivity comparison (CH4, C2_4, C5+)
    groups = ["CH4", "C2_4", "C5+"]
    sel_model = {}
    # compute outlet selectivity from F_model
    F_out = F_model[:, -1]
    prod_CH4 = max(F_out[species_keys["CH4"]], 0.0)
    prod_C2_4 = max(F_out[species_keys["C2_4"]], 0.0)
    prod_C5 = max(F_out[species_keys["C5plus"]], 0.0)
    total = prod_CH4 + prod_C2_4 + prod_C5
    if total > 0:
        sel_model = {"CH4": prod_CH4 / total, "C2_4": prod_C2_4 / total, "C5+": prod_C5 / total}
    else:
        sel_model = {g: 0.0 for g in groups}

    sel_exp = {"CH4": exp_record.get("sel_CH4", 0.0), "C2_4": exp_record.get("sel_C2_4", 0.0), "C5+": exp_record.get("sel_C5_plus", 0.0)}

    x = np.arange(len(groups))
    width = 0.35
    axes[1].bar(x - width/2, [sel_model[g] for g in groups], width, label="Model")
    axes[1].bar(x + width/2, [sel_exp[g] for g in groups], width, label="Exp")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(groups)
    axes[1].set_ylabel("Selectivity")
    axes[1].legend()

    # 3) Product distribution (optional: alkene/alkane fractions if present)
    # If record contains alkene_fraction and alkane_fraction arrays, plot them
    if exp_record.get("alkene_fraction") is not None and exp_record.get("alkane_fraction") is not None:
        af = exp_record["alkene_fraction"]
        alf = exp_record["alkane_fraction"]
        n = min(len(af), len(alf))
        x = np.arange(1, n+1)
        axes[2].plot(x, af[:n], label="Alkenes")
        axes[2].plot(x, alf[:n], label="Alkanes")
        axes[2].set_xlabel("Carbon number (group)")
        axes[2].set_ylabel("Fraction")
        axes[2].legend()
    else:
        axes[2].text(0.5, 0.5, "No product distribution data", ha='center', va='center')
        axes[2].set_axis_off()

    plt.tight_layout()
    if outpath:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, dpi=200)
        plt.close(fig)
    else:
        plt.show()
    return fig


# Diagnostics utilities

def compute_model_diagnostics(model, T, P, z, y):
    """Compute diagnostics (theta_CH2 and key rates) along axial positions.

    Returns a dict with arrays for: theta_CH2, r_growth (r4g), r_ch4 (r5), r_c2_4 (r6), r_c5p (r7)
    """
    import numpy as np

    n = y.shape[1]
    theta_CH2 = np.zeros(n)
    r_growth = np.zeros(n)
    r_ch4 = np.zeros(n)
    r_c2_4 = np.zeros(n)
    r_c5p = np.zeros(n)

    for i in range(n):
        Fi = y[:, i]
        # If model provides a solve_surface, use it for coverages
        if hasattr(model, "solve_surface"):
            cov = model.solve_surface(T, P, Fi)
            theta_CH2[i] = cov.get("theta_CH2", cov.get("theta_CH_2", 0.0) if cov is not None else 0.0)
        else:
            theta_CH2[i] = 0.0
        rates = model.rate(T, P, Fi)
        rates = np.atleast_1d(rates)
        # indices per BrubachModel: r4g index 4, r5 index 5, r6 index 6, r7 index 7
        if rates.size >= 8:
            r_growth[i] = rates[4]
            r_ch4[i] = rates[5]
            r_c2_4[i] = rates[6]
            r_c5p[i] = rates[7]
        else:
            # best-effort mapping for simpler models
            r_growth[i] = 0.0
            r_ch4[i] = rates[-1] if rates.size >= 1 else 0.0
            r_c2_4[i] = 0.0
            r_c5p[i] = 0.0

    return {
        "theta_CH2": theta_CH2,
        "r_growth": r_growth,
        "r_ch4": r_ch4,
        "r_c2_4": r_c2_4,
        "r_c5p": r_c5p,
    }


def plot_diagnostics(z, diag, outpath=None):
    """Plot theta_CH2 and rates vs reactor length."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.plot(z, diag["theta_CH2"], label=r"$\theta_{CH2}$")
    ax1.set_ylabel(r"$\theta_{CH2}$")
    ax1.grid(True)

    ax2.plot(z, diag["r_growth"], label="r_growth")
    ax2.plot(z, diag["r_ch4"], label="r_CH4")
    ax2.plot(z, diag["r_c2_4"], label="r_C2_4")
    ax2.plot(z, diag["r_c5p"], label="r_C5+")
    ax2.set_ylabel("Rate (mol m^-3 s^-1)")
    ax2.set_xlabel("Reactor length, m")
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    if outpath:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, dpi=200)
        plt.close(fig)
    else:
        plt.show()
    return fig
