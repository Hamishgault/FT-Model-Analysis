"""Plotting helpers for FT model comparisons.

These are simple wrappers around matplotlib to produce conversion and
selectivity plots. Kept small so they can be adapted or replaced.
"""
from typing import Dict, Iterable
import matplotlib.pyplot as plt
import numpy as np


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
