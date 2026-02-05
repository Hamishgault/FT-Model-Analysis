"""
Structural sensitivity analysis and operating-space scan for FT/RWGS reactor.

This script reuses the existing solve_until_equilibrium and extract_profiles
helpers without modifying the reactor equations.

Outputs:
  - sensitivity_table.csv
  - sensitivity_ranking.csv
  - operating_map.png
"""

import os
import sys
from copy import deepcopy
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT_ROOT)

from examples.compare_co2_feed_analysis import (
    SIM_CONFIG,
    PRODUCT_COMPONENTS,
    build_inlet_flow,
    solve_until_equilibrium,
    extract_profiles,
    solve_case_with_retry,
)

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "examples", "outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)


SENSITIVITY_PARAMS = [
    "k_rwgs",
    "k_c1",
    "k_c2_c4",
    "k_c5_c12",
    "k_c13_plus",
    "k_cracking",
    "k_light_cracking",
    "eta_ft",
    "ua_per_kg",
    "dp_dw",
]

METRICS = [
    "co2_conversion",
    "co_selectivity",
    "c5_c12_yield",
    "c13_plus_yield",
]


def safe_solve(sim_config: Dict[str, float], inlet_flow: Dict[str, float]):
    """Run solve_until_equilibrium safely; return profiles or None on failure."""
    try:
        return solve_until_equilibrium(sim_config, inlet_flow)
    except Exception:
        return None


def _equilibrium_index(W: np.ndarray, W_eq: float) -> int:
    if W.size == 0:
        return 0
    return int(np.argmin(np.abs(W - W_eq)))


def compute_metrics_from_profiles(
    profiles,
    inlet_flow: Dict[str, float],
) -> Dict[str, float]:
    """
    Compute metrics at equilibrium index using available mole fraction profiles.

    Notes
    -----
    Conversion/selectivity here are based on mole fractions because the helper
    profiles do not return total molar flow. This is suitable for screening.
    """
    (W, product_mole_frac, _mapped_mole_frac, co2_mole_frac, all_mole_frac, W_eq) = profiles

    idx = _equilibrium_index(W, W_eq)
    total_inlet = sum(inlet_flow.values())
    if total_inlet <= 0:
        return {m: np.nan for m in METRICS}

    co2_in = inlet_flow.get("CO2", 0.0) / total_inlet
    co2_out = co2_mole_frac[idx] / 100.0

    co2_conversion = np.nan
    if co2_in > 0:
        co2_conversion = max(0.0, min(1.0, (co2_in - co2_out) / co2_in))

    co_out = all_mole_frac.get("CO", np.array([np.nan]))[idx] / 100.0

    product_sum = 0.0
    for comp in PRODUCT_COMPONENTS:
        values = product_mole_frac.get(comp)
        if values is None:
            continue
        product_sum += values[idx] / 100.0

    denom = co_out + product_sum
    co_selectivity = np.nan if denom <= 0 else co_out / denom

    c5_c12_yield = np.nan
    if "C5_C12" in product_mole_frac:
        c5_c12_yield = product_mole_frac["C5_C12"][idx] / 100.0

    c13_plus_yield = np.nan
    if "C13_plus" in product_mole_frac:
        c13_plus_yield = product_mole_frac["C13_plus"][idx] / 100.0

    return {
        "co2_conversion": co2_conversion,
        "co_selectivity": co_selectivity,
        "c5_c12_yield": c5_c12_yield,
        "c13_plus_yield": c13_plus_yield,
    }


def run_single_case(sim_config: Dict[str, float], inlet_flow: Dict[str, float]) -> Tuple[Dict[str, float], Optional[Tuple]]:
    profiles = safe_solve(sim_config, inlet_flow)
    if profiles is None:
        return {m: np.nan for m in METRICS}, None
    return compute_metrics_from_profiles(profiles, inlet_flow), profiles


def sensitivity_screen(
    base_config: Dict[str, float],
    inlet_flow: Dict[str, float],
    perturb_fraction: float = 0.30,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run ± perturbations and compute normalized sensitivity indices.
    """
    records: List[Dict[str, float]] = []
    sens_records: List[Dict[str, float]] = []

    base_metrics, _base_profiles = run_single_case(base_config, inlet_flow)

    for param in SENSITIVITY_PARAMS:
        base_value = base_config.get(param)
        if base_value is None:
            continue

        for sign in (-1, 1):
            config = deepcopy(base_config)
            config[param] = base_value * (1.0 + sign * perturb_fraction)
            metrics, _profiles = run_single_case(config, inlet_flow)

            records.append(
                {
                    "parameter": param,
                    "perturbation": f"{sign * perturb_fraction:+.0%}",
                    "value": config[param],
                    "base_value": base_value,
                    **metrics,
                }
            )

        # Compute central-difference sensitivity index for each metric
        minus = next(
            r for r in records if r["parameter"] == param and r["perturbation"] == f"{-perturb_fraction:+.0%}"
        )
        plus = next(
            r for r in records if r["parameter"] == param and r["perturbation"] == f"{perturb_fraction:+.0%}"
        )

        for metric in METRICS:
            y_base = base_metrics.get(metric, np.nan)
            y_minus = minus.get(metric, np.nan)
            y_plus = plus.get(metric, np.nan)

            delta_y = y_plus - y_minus
            delta_p = 2.0 * perturb_fraction
            if not np.isfinite(y_base) or y_base == 0 or not np.isfinite(delta_y):
                s_val = np.nan
            else:
                s_val = (delta_y / y_base) / delta_p

            sens_records.append(
                {
                    "parameter": param,
                    "metric": metric,
                    "S": s_val,
                    "abs_S": np.nan if not np.isfinite(s_val) else abs(s_val),
                }
            )

    sensitivity_table = pd.DataFrame.from_records(records)
    sensitivity_ranking = pd.DataFrame.from_records(sens_records)

    if not sensitivity_ranking.empty:
        sensitivity_ranking = sensitivity_ranking.sort_values(
            ["metric", "abs_S"], ascending=[True, False]
        )

    return sensitivity_table, sensitivity_ranking


def operating_space_scan(
    base_config: Dict[str, float],
    temperature_range: Tuple[float, float] = (500.0, 600.0),
    ratio_range: Tuple[float, float] = (2.0, 6.0),
    n_temp: int = 11,
    n_ratio: int = 9,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    2D scan over temperature and H2/CO2 ratio.
    """
    temps = np.linspace(temperature_range[0], temperature_range[1], n_temp)
    ratios = np.linspace(ratio_range[0], ratio_range[1], n_ratio)

    conversion_map = np.full((n_ratio, n_temp), np.nan)
    c5_selectivity_map = np.full((n_ratio, n_temp), np.nan)

    for i, ratio in enumerate(ratios):
        co2_fraction = 1.0 / (1.0 + ratio)
        inlet_flow = build_inlet_flow(co2_fraction, total_flow=1.0)

        for j, temp in enumerate(temps):
            config = deepcopy(base_config)
            config["temperature"] = temp

            metrics, profiles = run_single_case(config, inlet_flow)
            conversion_map[i, j] = metrics["co2_conversion"]

            # C5–C12 selectivity based on product distribution
            if profiles is None:
                c5_selectivity_map[i, j] = np.nan
                continue

            (W, product_mole_frac, _mapped, _co2, _all, W_eq) = profiles
            idx = _equilibrium_index(W, W_eq)
            product_sum = 0.0
            for comp in PRODUCT_COMPONENTS:
                values = product_mole_frac.get(comp)
                if values is None:
                    continue
                product_sum += values[idx] / 100.0

            if product_sum <= 0:
                c5_selectivity_map[i, j] = np.nan
            else:
                c5_values = product_mole_frac.get("C5_C12")
                c5_selectivity_map[i, j] = (
                    np.nan if c5_values is None else (c5_values[idx] / 100.0) / product_sum
                )

    return temps, ratios, conversion_map, c5_selectivity_map


def plot_operating_map(
    temps: np.ndarray,
    ratios: np.ndarray,
    conversion_map: np.ndarray,
    c5_selectivity_map: np.ndarray,
    output_path: str,
) -> None:
    """Plot contour maps for conversion and C5–C12 selectivity."""
    T_grid, R_grid = np.meshgrid(temps, ratios)

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5), sharey=True)

    conv = axes[0].contourf(T_grid, R_grid, conversion_map, levels=20, cmap="viridis")
    axes[0].set_title("CO2 conversion")
    axes[0].set_xlabel("Temperature [K]")
    axes[0].set_ylabel("H2/CO2 ratio")
    fig.colorbar(conv, ax=axes[0])

    c5 = axes[1].contourf(T_grid, R_grid, c5_selectivity_map, levels=20, cmap="plasma")
    axes[1].set_title("C5–C12 selectivity")
    axes[1].set_xlabel("Temperature [K]")
    fig.colorbar(c5, ax=axes[1])

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)


def main():
    base_config = deepcopy(SIM_CONFIG)

    # Base inlet: H2/CO2 ratio = 3 (CO2 fraction = 0.25)
    base_ratio = 3.0
    base_inlet = build_inlet_flow(1.0 / (1.0 + base_ratio), total_flow=1.0)

    # One-time sanity extraction using the shared extract_profiles helper
    try:
        base_model, _used_config = solve_case_with_retry(base_config, base_inlet)
        if base_model is not None:
            _profiles = extract_profiles(base_model, base_config.get("equilibrium_tol", 1e-4))
    except Exception:
        pass

    sensitivity_table, sensitivity_ranking = sensitivity_screen(base_config, base_inlet)
    sensitivity_table.to_csv(os.path.join(OUTPUT_DIR, "sensitivity_table.csv"), index=False)
    sensitivity_ranking.to_csv(os.path.join(OUTPUT_DIR, "sensitivity_ranking.csv"), index=False)

    temps, ratios, conversion_map, c5_selectivity_map = operating_space_scan(base_config)
    plot_operating_map(
        temps,
        ratios,
        conversion_map,
        c5_selectivity_map,
        os.path.join(OUTPUT_DIR, "operating_map.png"),
    )

    print(f"[OK] Saved sensitivity_table.csv and sensitivity_ranking.csv to {OUTPUT_DIR}")
    print(f"[OK] Saved operating_map.png to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
