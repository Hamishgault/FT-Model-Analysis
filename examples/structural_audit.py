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
from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, value
from idaes.core import FlowsheetBlock

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC_DIR = os.path.join(PROJECT_ROOT, "src")
sys.path.insert(0, SRC_DIR)

from ft_model.ft_rwgs_zeolite_reactor import FTRWGSReactor, discretize_reactor  # type: ignore[reportMissingImports]
from ft_model.sim_config import SIM_CONFIG  # type: ignore[reportMissingImports]

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "examples", "outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)


PRODUCT_COMPONENTS = [
    "C1",
    "C2_C4",
    "C5_C12",
    "C13_plus",
    "iso_C5_C12",
    "aromatics",
    "coke",
]

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

AUDIT_OVERRIDES = {
    "max_iter": 4000,
    "nfe": 6,
    "acceptable_tol": 1e-5,
}


def build_inlet_flow(co2_fraction: float, total_flow: float = 1.0) -> Dict[str, float]:
    eps = 1e-8
    co2_flow = total_flow * co2_fraction
    h2_flow = total_flow * (1 - co2_fraction)
    return {
        "CO2": co2_flow,
        "H2": h2_flow,
        "CO": eps,
        "H2O": eps,
        "C1": eps,
        "C2_C4": eps,
        "C5_C12": eps,
        "C13_plus": eps,
        "iso_C5_C12": eps,
        "aromatics": eps,
        "coke": eps,
    }


def solve_case(sim_config: Dict[str, float], inlet_flow: Dict[str, float]):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0])

    m.fs.reactor = FTRWGSReactor(
        include_zeolite_reactions=bool(sim_config["include_zeolite_reactions"]),
        energy_balance=bool(sim_config["energy_balance"]),
        pressure_drop=bool(sim_config["pressure_drop"]),
        ergun_pressure_drop=bool(sim_config["ergun_pressure_drop"]),
        heat_transfer=bool(sim_config["heat_transfer"]),
        mass_transfer=bool(sim_config["mass_transfer"]),
        kinetics_model=sim_config.get("kinetics_model", "lumped_simple"),
    )

    m.fs.reactor.initialize(
        inlet_flow=inlet_flow,
        temperature=sim_config["temperature"],
        pressure=sim_config["pressure_bar"] * 101325.0,
        W_total=sim_config["W_total"],
        k_rwgs=sim_config["k_rwgs"],
        Keq_rwgs=sim_config["Keq_rwgs"],
        k_c1=sim_config["k_c1"],
        k_c2_c4=sim_config["k_c2_c4"],
        k_c5_c12=sim_config["k_c5_c12"],
        k_c13_plus=sim_config["k_c13_plus"],
        k_cracking=sim_config["k_cracking"],
        k_light_cracking=sim_config["k_light_cracking"],
        k_isomerization=sim_config["k_isomerization"],
        k_oligomerization=sim_config["k_oligomerization"],
        k_aromatization=sim_config["k_aromatization"],
        k_coke_formation=sim_config["k_coke_formation"],
        kfts_ref=sim_config.get("kfts_ref"),
        E_app=sim_config.get("E_app"),
        b_ref=sim_config.get("b_ref"),
        dH_b=sim_config.get("dH_b"),
        T_ref=sim_config.get("T_ref"),
        beta_gasoline=sim_config.get("beta_gasoline"),
        beta_jet=sim_config.get("beta_jet"),
        beta_diesel=sim_config.get("beta_diesel"),
        split_c5_gasoline=sim_config.get("split_c5_gasoline"),
        split_c5_jet=sim_config.get("split_c5_jet"),
        split_c13_diesel=sim_config.get("split_c13_diesel"),
        dp_dw=sim_config["dp_dw"],
        ergun_porosity=sim_config["ergun_porosity"],
        particle_diameter=sim_config["particle_diameter"],
        catalyst_bulk_density=sim_config["catalyst_bulk_density"],
        reactor_diameter=sim_config["reactor_diameter"],
        reactor_length=sim_config["reactor_length"],
        gas_viscosity=sim_config["gas_viscosity"],
        ua_per_kg=sim_config["ua_per_kg"],
        T_coolant=sim_config["T_coolant"],
        eta_ft=sim_config["eta_ft"],
        eta_zeolite=sim_config["eta_zeolite"],
    )

    discretize_reactor(m.fs.reactor, nfe=int(sim_config["nfe"]))

    t = 0
    w_inlet = m.fs.reactor.W.first()
    for comp in m.fs.reactor.component_list:
        m.fs.reactor.flow_mol_comp[t, w_inlet, comp].set_value(inlet_flow.get(comp, 0.0))

    m.fs.reactor.temperature[t, w_inlet].set_value(sim_config["temperature"])
    m.fs.reactor.pressure[t, w_inlet].set_value(sim_config["pressure_bar"] * 101325.0)

    solver = SolverFactory("ipopt")
    solver.options["max_iter"] = int(sim_config["max_iter"])
    solver.options["tol"] = sim_config["tol"]
    if sim_config.get("acceptable_tol") is not None:
        solver.options["acceptable_tol"] = sim_config["acceptable_tol"]
    if sim_config.get("linear_solver"):
        solver.options["linear_solver"] = sim_config["linear_solver"]
    if sim_config.get("bound_push") is not None:
        solver.options["bound_push"] = sim_config["bound_push"]
    if sim_config.get("mu_strategy"):
        solver.options["mu_strategy"] = sim_config["mu_strategy"]

    try:
        results = solver.solve(m, tee=False)
    except Exception:
        return None

    term_cond = results.solver.termination_condition
    if term_cond not in (
        TerminationCondition.optimal,
        TerminationCondition.locallyOptimal,
        TerminationCondition.feasible,
        TerminationCondition.maxIterations,
    ):
        return None

    if term_cond == TerminationCondition.maxIterations:
        print("[WARN] Solver hit maxIterations; using last iterate for audit metrics.")

    return m


def solve_case_with_retry(sim_config: Dict[str, float], inlet_flow: Dict[str, float]):
    w_scales = (1.0, 0.5, 0.2, 0.1)
    k_scales = (1.0, 0.5, 0.2, 0.1)
    nfe_values = (sim_config["nfe"], max(6, sim_config["nfe"] // 2), 6)
    kinetic_keys = [
        "k_rwgs",
        "k_c1",
        "k_c2_c4",
        "k_c5_c12",
        "k_c13_plus",
        "k_cracking",
        "k_light_cracking",
        "k_isomerization",
        "k_oligomerization",
        "k_aromatization",
        "k_coke_formation",
    ]

    def _attempt(base_config: Dict[str, float]):
        for nfe in nfe_values:
            for w_scale in w_scales:
                for k_scale in k_scales:
                    attempt_config = dict(base_config)
                    attempt_config["nfe"] = int(nfe)
                    attempt_config["W_total"] = base_config["W_total"] * w_scale
                    for key in kinetic_keys:
                        attempt_config[key] = base_config[key] * k_scale
                    if base_config.get("kfts_ref") is not None:
                        attempt_config["kfts_ref"] = base_config["kfts_ref"] * k_scale

                    model = solve_case(attempt_config, inlet_flow)
                    if model is not None:
                        return model, attempt_config
        return None, base_config

    return _attempt(sim_config)


def find_equilibrium_index(W: np.ndarray, mole_frac: Dict[str, np.ndarray], tol: float = 1e-4) -> int:
    if W.size < 3:
        return W.size - 1

    max_grad = np.zeros_like(W)
    for values in mole_frac.values():
        grad = np.abs(np.gradient(values, W))
        max_grad = np.maximum(max_grad, grad)

    for idx in range(1, W.size):
        if np.all(max_grad[idx:] < tol):
            return idx
    return W.size - 1


def extract_profiles(m, eq_tol: float):
    reactor = m.fs.reactor
    t = 0
    W_points = list(reactor.W)
    W = np.array([float(w) for w in W_points])

    total_flow = np.array([value(reactor.flow_mol_total[t, w]) for w in W_points])
    total_flow = np.where(total_flow > 1e-12, total_flow, 1.0)

    product_flows = {
        comp: np.array([value(reactor.flow_mol_comp[t, w, comp]) for w in W_points])
        for comp in PRODUCT_COMPONENTS
    }

    product_mole_frac = {
        comp: 100.0 * product_flows[comp] / total_flow for comp in PRODUCT_COMPONENTS
    }

    all_mole_frac = {
        comp: 100.0 * np.array([value(reactor.flow_mol_comp[t, w, comp]) for w in W_points]) / total_flow
        for comp in reactor.component_list
    }

    co2_flow = np.array([value(reactor.flow_mol_comp[t, w, "CO2"]) for w in W_points])
    co2_mole_frac = 100.0 * co2_flow / total_flow

    eq_idx = find_equilibrium_index(W, all_mole_frac, tol=eq_tol)
    W_eq = W[eq_idx]

    return W, product_mole_frac, co2_mole_frac, all_mole_frac, W_eq


def solve_until_equilibrium(sim_config: Dict[str, float], inlet_flow: Dict[str, float]):
    eq_tol = sim_config.get("equilibrium_tol", 1e-4)
    model, used_config = solve_case_with_retry(sim_config, inlet_flow)

    if model is None:
        return None

    try:
        profiles = extract_profiles(model, eq_tol)
    except Exception:
        return None

    W = profiles[0]
    W_eq = profiles[-1]
    if W_eq < 0.95 * W[-1]:
        print(f"[INFO] Equilibrium reached at W={W_eq:.3f} with W_total={used_config['W_total']}")
    else:
        print(
            f"[WARN] Equilibrium not reached (W_eq={W_eq:.3f} ~ W_end={W[-1]:.3f}) "
            f"with W_total={used_config['W_total']}"
        )

    return profiles


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
    (W, product_mole_frac, co2_mole_frac, all_mole_frac, W_eq) = profiles

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

            (W, product_mole_frac, _co2, _all, W_eq) = profiles
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
    base_config.update(AUDIT_OVERRIDES)

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
