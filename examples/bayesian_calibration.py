"""
Bayesian calibration for FT/RWGS reactor parameters using emcee + ArviZ.

Expected experimental data CSV columns:
- temperature (K)
- pressure_bar (bar)
- h2_co2_ratio
- co2_conversion
- co_selectivity
- c5_c12_yield
- c13_plus_yield
Optional uncertainty columns (one per metric, 1-sigma):
- sigma_co2_conversion
- sigma_co_selectivity
- sigma_c5_c12_yield
- sigma_c13_plus_yield

Notes:
- This script reuses solve_until_equilibrium from compare_co2_feed_analysis.
- It does NOT change reactor equations.
"""

import argparse
import os
import sys
import time
from multiprocessing import Pool, cpu_count
from copy import deepcopy
from typing import Dict, List

import numpy as np
import pandas as pd
import emcee
import arviz as az

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT_ROOT)

from src.ft_model.ft_rwgs_zeolite_reactor import FTRWGSReactor, discretize_reactor


METRICS = [
    "co2_conversion",
    "co_selectivity",
    "c5_c12_yield",
    "c13_plus_yield",
]

PARAM_NAMES = [
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

BASE_CONFIG: Dict[str, float] = {
    # Model toggles
    "include_zeolite_reactions": True,
    "energy_balance": False,
    "pressure_drop": False,
    "ergun_pressure_drop": False,
    "heat_transfer": False,
    "mass_transfer": False,
    "kinetics_model": "marvast_2005",

    # Operating conditions
    "temperature": 523.15,
    "pressure_bar": 23.0,
    "W_total": 5.0,
    "nfe": 8,

    # Kinetics (base case)
    "k_rwgs": 1e-4,
    "Keq_rwgs": 0.8,
    "k_c1": 1e-4,
    "k_c2_c4": 1e-4,
    "k_c5_c12": 1e-4,
    "k_c13_plus": 1e-4,
    "k_cracking": 1e-5,
    "k_light_cracking": 1e-5,
    "k_isomerization": 1e-5,
    "k_oligomerization": 8e-6,
    "k_aromatization": 5e-6,
    "k_coke_formation": 1e-6,
    "kfts_ref": 6.4e-4,
    "E_app": 23000.0,
    "b_ref": 1.6e-2,
    "dH_b": -28500.0,
    "T_ref": 543.0,
    "dp_dw": 1.0e3,
    "ergun_porosity": 0.40,
    "particle_diameter": 5.0e-3,
    "catalyst_bulk_density": 1000.0,
    "reactor_diameter": 1.0,
    "reactor_length": 1.0,
    "gas_viscosity": 1.0e-5,
    "ua_per_kg": 0.01,
    "T_coolant": 500.0,
    "eta_ft": 0.9,
    "eta_zeolite": 0.8,

    # Solver
    "max_iter": 600,
    "tol": 1e-6,
    "linear_solver": "mumps",
    "staged_solve": False,
}

FAST_OVERRIDES: Dict[str, float] = {
    "nfe": 6,
    "W_total": 2.0,
    "max_iter": 300,
    "staged_solve": False,
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


def safe_solve(sim_config: Dict[str, float], inlet_flow: Dict[str, float]):
    """Solve a single model instance and return profiles or None on failure."""
    try:
        model = solve_case(sim_config, inlet_flow)
    except Exception:
        return None
    if model is None:
        return None

    return extract_profiles(model)


def _equilibrium_index(W: np.ndarray, W_eq: float) -> int:
    if W.size == 0:
        return 0
    return int(np.argmin(np.abs(W - W_eq)))


def compute_metrics_from_profiles(profiles, inlet_flow: Dict[str, float]) -> Dict[str, float]:
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
    for comp, values in product_mole_frac.items():
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


def run_model_metrics(sim_config: Dict[str, float], row: pd.Series) -> Dict[str, float]:
    inlet_flow = build_inlet_flow(1.0 / (1.0 + row["h2_co2_ratio"]), total_flow=1.0)
    config = deepcopy(sim_config)
    config["temperature"] = float(row["temperature"])
    config["pressure_bar"] = float(row["pressure_bar"])

    profiles = safe_solve(config, inlet_flow)
    if profiles is None:
        return {m: np.nan for m in METRICS}

    return compute_metrics_from_profiles(profiles, inlet_flow)


def log_prior(theta: np.ndarray, base_config: Dict[str, float], param_names: List[str]) -> float:
    # Log-normal priors for positive rate/transport constants
    logp = 0.0
    for i, name in enumerate(param_names):
        val = theta[i]
        base_val = base_config.get(name)
        if base_val is None:
            return -np.inf

        if name == "eta_ft":
            if val <= 0.0 or val >= 1.5:
                return -np.inf
            mu, sigma = 0.9, 0.2
            logp += -0.5 * ((val - mu) / sigma) ** 2
            continue

        if val <= 0:
            return -np.inf

        # log-normal around base value, 30% 1-sigma in log space
        mu = np.log(base_val)
        sigma = 0.3
        logp += -0.5 * ((np.log(val) - mu) / sigma) ** 2

    return logp


def log_likelihood(
    theta: np.ndarray,
    base_config: Dict[str, float],
    data: pd.DataFrame,
    param_names: List[str],
) -> float:
    config = deepcopy(base_config)
    for i, name in enumerate(param_names):
        config[name] = float(theta[i])

    ll = 0.0
    for _, row in data.iterrows():
        preds = run_model_metrics(config, row)

        for metric in METRICS:
            y = float(row[metric])
            yhat = preds.get(metric, np.nan)
            if not np.isfinite(yhat):
                return -np.inf

            sigma_col = f"sigma_{metric}"
            sigma = float(row[sigma_col]) if sigma_col in row and pd.notna(row[sigma_col]) else 0.02
            if sigma <= 0:
                return -np.inf

            ll += -0.5 * ((yhat - y) / sigma) ** 2 - np.log(sigma * np.sqrt(2 * np.pi))

    return ll


def log_posterior(
    theta: np.ndarray,
    base_config: Dict[str, float],
    data: pd.DataFrame,
    param_names: List[str],
) -> float:
    lp = log_prior(theta, base_config, param_names)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta, base_config, data, param_names)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


def run_mcmc(
    data_csv: str,
    nwalkers: int = 24,
    nsteps: int = 200,
    burn_in: int = 50,
    nworkers: int | None = None,
    report_every: int = 10,
    param_names: List[str] | None = None,
    fast_mode: bool = True,
) -> None:
    data = pd.read_csv(data_csv)
    base_config = deepcopy(BASE_CONFIG)
    if fast_mode:
        base_config.update(FAST_OVERRIDES)

    if param_names is None:
        param_names = PARAM_NAMES

    ndim = len(param_names)
    base_theta = np.array([base_config[name] for name in param_names], dtype=float)

    # Initialize walkers around base values
    rng = np.random.default_rng(123)
    pos = base_theta * (1.0 + 0.05 * rng.standard_normal((nwalkers, ndim)))

    if nworkers is None:
        nworkers = max(1, cpu_count() - 1)

    with Pool(processes=nworkers) as pool:
        sampler = emcee.EnsembleSampler(
            nwalkers,
            ndim,
            log_posterior,
            args=(base_config, data, param_names),
            pool=pool,
        )

        start = time.perf_counter()
        for step, _state in enumerate(sampler.sample(pos, iterations=nsteps, progress=False), start=1):
            if step % report_every == 0 or step == nsteps:
                elapsed = time.perf_counter() - start
                rate = elapsed / step
                remaining = rate * (nsteps - step)
                print(
                    f"[MCMC] Step {step}/{nsteps} | "
                    f"elapsed={elapsed:.1f}s | ETA={remaining:.1f}s"
                )

    samples = sampler.get_chain(discard=burn_in, flat=True)

    idata = az.from_emcee(sampler, var_names=param_names)
    az.to_netcdf(idata, os.path.join(PROJECT_ROOT, "examples", "outputs", "posterior.nc"))

    summary = az.summary(idata)
    summary.to_csv(os.path.join(PROJECT_ROOT, "examples", "outputs", "posterior_summary.csv"))

    print("[OK] Saved posterior.nc and posterior_summary.csv to examples/outputs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Bayesian calibration with emcee")
    parser.add_argument(
        "--data",
        default=os.path.join(PROJECT_ROOT, "examples", "data", "experimental_results.csv"),
        help="Path to experimental results CSV",
    )
    parser.add_argument(
        "--params",
        default=",".join(PARAM_NAMES),
        help="Comma-separated parameter names to calibrate",
    )
    parser.add_argument("--nwalkers", type=int, default=24)
    parser.add_argument("--nsteps", type=int, default=200)
    parser.add_argument("--burn-in", type=int, default=50)
    parser.add_argument("--nworkers", type=int, default=None)
    parser.add_argument("--report-every", type=int, default=10)
    parser.add_argument("--fast", action="store_true", default=True)
    parser.add_argument("--no-fast", dest="fast", action="store_false")
    args = parser.parse_args()

    if not os.path.exists(args.data):
        raise FileNotFoundError(
            f"Experimental data not found: {args.data}. "
            "Create the CSV with required columns before running."
        )

    param_names = [p.strip() for p in args.params.split(",") if p.strip()]
    run_mcmc(
        args.data,
        nwalkers=args.nwalkers,
        nsteps=args.nsteps,
        burn_in=args.burn_in,
        nworkers=args.nworkers,
        report_every=args.report_every,
        param_names=param_names,
        fast_mode=args.fast,
    )


def solve_case(sim_config: Dict[str, float], inlet_flow: Dict[str, float]):
    from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition
    from idaes.core import FlowsheetBlock

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
    W_inlet = m.fs.reactor.W.first()
    for c in m.fs.reactor.component_list:
        m.fs.reactor.flow_mol_comp[t, W_inlet, c].set_value(inlet_flow.get(c, 0.0))

    m.fs.reactor.temperature[t, W_inlet].set_value(sim_config["temperature"])
    m.fs.reactor.pressure[t, W_inlet].set_value(sim_config["pressure_bar"] * 101325.0)

    solver = SolverFactory("ipopt")
    solver.options["max_iter"] = int(sim_config["max_iter"])
    solver.options["tol"] = sim_config["tol"]
    if sim_config.get("linear_solver"):
        solver.options["linear_solver"] = sim_config["linear_solver"]

    try:
        results = solver.solve(m, tee=False, load_solutions=False)
    except Exception:
        return None

    if results.solver.termination_condition != TerminationCondition.optimal:
        return None

    try:
        m.solutions.load_from(results)
    except ValueError:
        return None

    return m


def extract_profiles(m) -> tuple:
    reactor = m.fs.reactor
    t = 0
    W_points = list(reactor.W)
    W = np.array([float(w) for w in W_points])

    total_flow = np.array([reactor.flow_mol_total[t, w].value for w in W_points])
    total_flow = np.where(total_flow > 1e-12, total_flow, 1.0)

    product_components = ["C1", "C2_C4", "C5_C12", "C13_plus", "iso_C5_C12", "aromatics", "coke"]
    product_flows = {
        comp: np.array([reactor.flow_mol_comp[t, w, comp].value for w in W_points])
        for comp in product_components
    }
    product_mole_frac = {comp: 100.0 * product_flows[comp] / total_flow for comp in product_components}

    all_mole_frac = {
        comp: 100.0 * np.array([reactor.flow_mol_comp[t, w, comp].value for w in W_points]) / total_flow
        for comp in reactor.component_list
    }

    co2_flow = np.array([reactor.flow_mol_comp[t, w, "CO2"].value for w in W_points])
    co2_mole_frac = 100.0 * co2_flow / total_flow

    W_eq = W[-1]
    return W, product_mole_frac, {}, co2_mole_frac, all_mole_frac, W_eq
