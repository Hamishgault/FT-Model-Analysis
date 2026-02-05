import logging
import os
import sys
from typing import Dict, Tuple, List

import numpy as np
import matplotlib.pyplot as plt
from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, value
from idaes.core import FlowsheetBlock

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "src"))

from ft_model.ft_rwgs_zeolite_reactor import FTRWGSReactor, discretize_reactor  # type: ignore[reportMissingImports]

logging.getLogger("pyomo.repn.plugins.nl_writer").setLevel(logging.ERROR)

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "examples", "outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)


SIM_CONFIG: Dict[str, float] = {
    # Model toggles
    'include_zeolite_reactions': True,
    'energy_balance': True,
    'pressure_drop': True,
    'ergun_pressure_drop': True,
    'heat_transfer': True,
    'mass_transfer': True,
    'kinetics_model': 'rwgs_2017',

    # Operating conditions
    'temperature': 523.15,   # K
    'pressure_bar': 23.0,    # bar
    'W_total': 100.0,          # kg catalyst
    'nfe': 200,               # spatial elements

    # Kinetics (base case)
    'k_rwgs': 0.001,
    'Keq_rwgs': 0.8,
    'k_c1': 0.0002,
    'k_c2_c4': 0.0001,
    'k_c5_c12': 0.00005,
    'k_c13_plus': 0.00002,
    'k_cracking': 0.00002,
    'k_light_cracking': 0.00001,
    'k_isomerization': 0.00001,
    'k_oligomerization': 0.000008,
    'k_aromatization': 0.000005,
    'k_coke_formation': 0.000001,
    'kfts_ref': 6.4e-4,
    'E_app': 23000.0,
    'b_ref': 1.6e-2,
    'dH_b': -28500.0,
    'T_ref': 543.0,
    'dp_dw': 1.0e3,
    'ergun_porosity': 0.40,
    'particle_diameter': 5.0e-3,
    'catalyst_bulk_density': 1000.0,
    'reactor_diameter': 1.0,
    'reactor_length': 1.0,
    'gas_viscosity': 1.0e-5,
    'ua_per_kg': 0.01,
    'T_coolant': 500.0,
    'eta_ft': 0.9,
    'eta_zeolite': 0.8,

    # Selectivity/cracking factors (defaults)
    'beta_gasoline': 0.7,
    'beta_jet': 0.8,
    'beta_diesel': 0.6,
    'split_c5_gasoline': 0.5,
    'split_c5_jet': 0.5,
    'split_c13_diesel': 0.7,

    # Solver
    'max_iter': 2000,
    'tol': 1e-6,
    'acceptable_tol': 1e-5,
    'linear_solver': 'mumps',
    'bound_push': 1e-8,
    'mu_strategy': 'adaptive',
    'staged_solve': True,

    # Equilibrium detection
    'equilibrium_tol': 5e-4,
}

SIM_CONFIG.update({
    'energy_balance': False,
    'pressure_drop': False,
    'heat_transfer': False,
    'mass_transfer': False,
    'nfe': 8,
    'W_total': 5.0,
    'k_rwgs': 1e-4,
    'k_c1': 1e-4,
    'k_c2_c4': 1e-4,
    'k_c5_c12': 1e-4,
    'k_c13_plus': 1e-4,
    'k_cracking': 1e-5,
})

PRODUCT_COMPONENTS = [
    'C1',
    'C2_C4',
    'C5_C12',
    'C13_plus',
    'iso_C5_C12',
    'aromatics',
    'coke',
]

MAPPED_COMPONENTS = [
    'gasoline',
    'jet',
    'diesel',
]


def build_inlet_flow(co2_fraction: float, total_flow: float = 1.0) -> Dict[str, float]:
    eps = 1e-8
    co2_flow = total_flow * co2_fraction
    h2_flow = total_flow * (1 - co2_fraction)
    return {
        'CO2': co2_flow,
        'H2': h2_flow,
        'CO': eps,
        'H2O': eps,
        'C1': eps,
        'C2_C4': eps,
        'C5_C12': eps,
        'C13_plus': eps,
        'iso_C5_C12': eps,
        'aromatics': eps,
        'coke': eps,
    }


def solve_case(sim_config: Dict[str, float], inlet_flow: Dict[str, float]):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0])

    m.fs.reactor = FTRWGSReactor(
        include_zeolite_reactions=bool(sim_config['include_zeolite_reactions']),
        energy_balance=bool(sim_config['energy_balance']),
        pressure_drop=bool(sim_config['pressure_drop']),
        ergun_pressure_drop=bool(sim_config['ergun_pressure_drop']),
        heat_transfer=bool(sim_config['heat_transfer']),
        mass_transfer=bool(sim_config['mass_transfer']),
        kinetics_model=sim_config.get('kinetics_model', 'lumped_simple'),
    )

    m.fs.reactor.initialize(
        inlet_flow=inlet_flow,
        temperature=sim_config['temperature'],
        pressure=sim_config['pressure_bar'] * 101325.0,
        W_total=sim_config['W_total'],
        k_rwgs=sim_config['k_rwgs'],
        Keq_rwgs=sim_config['Keq_rwgs'],
        k_c1=sim_config['k_c1'],
        k_c2_c4=sim_config['k_c2_c4'],
        k_c5_c12=sim_config['k_c5_c12'],
        k_c13_plus=sim_config['k_c13_plus'],
        k_cracking=sim_config['k_cracking'],
        k_light_cracking=sim_config['k_light_cracking'],
        k_isomerization=sim_config['k_isomerization'],
        k_oligomerization=sim_config['k_oligomerization'],
        k_aromatization=sim_config['k_aromatization'],
        k_coke_formation=sim_config['k_coke_formation'],
        kfts_ref=sim_config.get('kfts_ref'),
        E_app=sim_config.get('E_app'),
        b_ref=sim_config.get('b_ref'),
        dH_b=sim_config.get('dH_b'),
        T_ref=sim_config.get('T_ref'),
        beta_gasoline=sim_config.get('beta_gasoline'),
        beta_jet=sim_config.get('beta_jet'),
        beta_diesel=sim_config.get('beta_diesel'),
        split_c5_gasoline=sim_config.get('split_c5_gasoline'),
        split_c5_jet=sim_config.get('split_c5_jet'),
        split_c13_diesel=sim_config.get('split_c13_diesel'),
        dp_dw=sim_config['dp_dw'],
        ergun_porosity=sim_config['ergun_porosity'],
        particle_diameter=sim_config['particle_diameter'],
        catalyst_bulk_density=sim_config['catalyst_bulk_density'],
        reactor_diameter=sim_config['reactor_diameter'],
        reactor_length=sim_config['reactor_length'],
        gas_viscosity=sim_config['gas_viscosity'],
        ua_per_kg=sim_config['ua_per_kg'],
        T_coolant=sim_config['T_coolant'],
        eta_ft=sim_config['eta_ft'],
        eta_zeolite=sim_config['eta_zeolite'],
    )

    discretize_reactor(m.fs.reactor, nfe=int(sim_config['nfe']))

    t = 0
    W_inlet = m.fs.reactor.W.first()
    for c in m.fs.reactor.component_list:
        m.fs.reactor.flow_mol_comp[t, W_inlet, c].set_value(inlet_flow.get(c, 0.0))

    m.fs.reactor.temperature[t, W_inlet].set_value(sim_config['temperature'])
    m.fs.reactor.pressure[t, W_inlet].set_value(sim_config['pressure_bar'] * 101325.0)

    solver = SolverFactory('ipopt')
    solver.options['max_iter'] = int(sim_config['max_iter'])
    solver.options['tol'] = sim_config['tol']
    if sim_config.get('acceptable_tol') is not None:
        solver.options['acceptable_tol'] = sim_config['acceptable_tol']
    if sim_config.get('linear_solver'):
        solver.options['linear_solver'] = sim_config['linear_solver']
    if sim_config.get('bound_push') is not None:
        solver.options['bound_push'] = sim_config['bound_push']
    if sim_config.get('mu_strategy'):
        solver.options['mu_strategy'] = sim_config['mu_strategy']

    if sim_config.get('staged_solve', True):
        original_eta_ft = m.fs.reactor.eta_ft.value
        original_eta_zeo = m.fs.reactor.eta_zeolite.value
        original_ua = m.fs.reactor.ua_per_kg.value
        original_rate_multiplier = m.fs.reactor.rate_multiplier.value
        original_pressure_drop_multiplier = m.fs.reactor.pressure_drop_multiplier.value

        if m.fs.reactor.config.mass_transfer:
            m.fs.reactor.eta_ft.set_value(1.0)
            m.fs.reactor.eta_zeolite.set_value(1.0)

        if m.fs.reactor.config.heat_transfer:
            m.fs.reactor.ua_per_kg.set_value(0.0)

        if m.fs.reactor.config.pressure_drop:
            m.fs.reactor.pressure_drop_multiplier.set_value(0.0)

        for ramp_value in (0.0, 0.01, 0.03, 0.1, 0.3, 0.6, 1.0):
            m.fs.reactor.rate_multiplier.set_value(ramp_value)
            solver.options['max_iter'] = 400
            stage_results = solver.solve(m, tee=False)
            if stage_results.solver.termination_condition != TerminationCondition.optimal:
                break

        if m.fs.reactor.config.pressure_drop:
            for pd_value in (0.0, 0.1, 0.3, 0.6, 1.0):
                m.fs.reactor.pressure_drop_multiplier.set_value(pd_value)
                solver.options['max_iter'] = 400
                stage_results = solver.solve(m, tee=False)
                if stage_results.solver.termination_condition != TerminationCondition.optimal:
                    break

        m.fs.reactor.eta_ft.set_value(original_eta_ft)
        m.fs.reactor.eta_zeolite.set_value(original_eta_zeo)
        m.fs.reactor.ua_per_kg.set_value(original_ua)
        m.fs.reactor.rate_multiplier.set_value(original_rate_multiplier)
        m.fs.reactor.pressure_drop_multiplier.set_value(original_pressure_drop_multiplier)
        solver.options['max_iter'] = int(sim_config['max_iter'])

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


def solve_case_with_retry(sim_config: Dict[str, float], inlet_flow: Dict[str, float]):
    w_scales = (1.0, 0.5, 0.2, 0.1, 0.05)
    k_scales = (1.0, 0.5, 0.2, 0.1, 0.03, 0.01)
    nfe_values = (sim_config['nfe'], max(6, sim_config['nfe'] // 2), 6)
    kinetic_keys = [
        'k_rwgs', 'k_c1', 'k_c2_c4', 'k_c5_c12', 'k_c13_plus',
        'k_cracking', 'k_light_cracking', 'k_isomerization',
        'k_oligomerization', 'k_aromatization', 'k_coke_formation',
    ]

    ua_scales = (1.0, 0.5, 0.2)
    dp_scales = (1.0, 0.5, 0.2)

    def _attempt(base_config: Dict[str, float]):
        for nfe in nfe_values:
            for w_scale in w_scales:
                for k_scale in k_scales:
                    for ua_scale in ua_scales:
                        for dp_scale in dp_scales:
                            attempt_config = dict(base_config)
                            attempt_config['nfe'] = int(nfe)
                            attempt_config['W_total'] = base_config['W_total'] * w_scale
                            attempt_config['ua_per_kg'] = base_config['ua_per_kg'] * ua_scale
                            attempt_config['dp_dw'] = base_config['dp_dw'] * dp_scale
                            for key in kinetic_keys:
                                attempt_config[key] = base_config[key] * k_scale
                            if 'kfts_ref' in base_config and base_config.get('kfts_ref') is not None:
                                attempt_config['kfts_ref'] = base_config['kfts_ref'] * k_scale

                            model = solve_case(attempt_config, inlet_flow)
                            if model is not None:
                                return model, attempt_config
        return None, base_config

    model, config_used = _attempt(sim_config)
    if model is not None:
        return model, config_used

    softened = dict(sim_config)
    softened['pressure_drop'] = False
    model, config_used = _attempt(softened)
    if model is not None:
        return model, config_used

    softened['heat_transfer'] = False
    model, config_used = _attempt(softened)
    if model is not None:
        return model, config_used

    softened['include_zeolite_reactions'] = False
    model, config_used = _attempt(softened)
    if model is not None:
        return model, config_used

    return None, config_used


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


def extract_profiles(m, eq_tol: float) -> Tuple[np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], np.ndarray, Dict[str, np.ndarray], float]:
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

    mapped_flows = {
        'gasoline': np.array([value(reactor.flow_gasoline[t, w]) for w in W_points]),
        'jet': np.array([value(reactor.flow_jet[t, w]) for w in W_points]),
        'diesel': np.array([value(reactor.flow_diesel[t, w]) for w in W_points]),
    }
    mapped_mole_frac = {
        comp: 100.0 * mapped_flows[comp] / total_flow for comp in MAPPED_COMPONENTS
    }

    co2_flow = np.array([value(reactor.flow_mol_comp[t, w, 'CO2']) for w in W_points])
    co2_mole_frac = 100.0 * co2_flow / total_flow

    eq_idx = find_equilibrium_index(W, all_mole_frac, tol=eq_tol)
    W_eq = W[eq_idx]

    return W, product_mole_frac, mapped_mole_frac, co2_mole_frac, all_mole_frac, W_eq


def solve_until_equilibrium(sim_config: Dict[str, float], inlet_flow: Dict[str, float]):
    eq_tol = sim_config.get('equilibrium_tol', 1e-4)
    try:
        model, used_config = solve_case_with_retry(sim_config, inlet_flow)
    except Exception:
        return None

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


def safe_forward(sim_config: Dict[str, float], inlet_flow: Dict[str, float]):
    try:
        return solve_until_equilibrium(sim_config, inlet_flow)
    except Exception:
        return None


def plot_results(
    results_by_case: Dict[str, Tuple[np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], np.ndarray, Dict[str, np.ndarray], float]]
):
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 10), sharex=True)

    for ax, (case_label, (W, product_mole_frac, _mapped_mole_frac, _co2, _all_mole_frac, W_eq)) in zip(axes, results_by_case.items()):
        for comp in PRODUCT_COMPONENTS:
            ax.plot(W, product_mole_frac[comp], label=comp)
        ax.set_title(f"Product distribution (mole %) vs reactor length ({case_label})")
        ax.set_ylabel("Mole fraction [%]")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=9)
        ax.axvline(W_eq, color="black", linestyle="--", alpha=0.5, label="Equilibrium")

    axes[-1].set_xlabel("Normalized catalyst weight W")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "products_mole_percent.png"), dpi=200)

    fig2, axes2 = plt.subplots(nrows=2, ncols=1, figsize=(10, 10), sharex=True)
    for ax, (case_label, (W, _product_mole_frac, mapped_mole_frac, _co2, _all_mole_frac, W_eq)) in zip(axes2, results_by_case.items()):
        for comp in MAPPED_COMPONENTS:
            ax.plot(W, mapped_mole_frac[comp], label=comp)
        ax.set_title(f"Mapped gasoline/jet/diesel (mole %) vs reactor length ({case_label})")
        ax.set_ylabel("Mole fraction [%]")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=9)
        ax.axvline(W_eq, color="black", linestyle="--", alpha=0.5, label="Equilibrium")

    axes2[-1].set_xlabel("Normalized catalyst weight W")
    fig2.tight_layout()
    fig2.savefig(os.path.join(OUTPUT_DIR, "mapped_fractions_mole_percent.png"), dpi=200)

    fig3, ax3 = plt.subplots(figsize=(10, 5))
    for case_label, (W, _product_mole_frac, _mapped_mole_frac, co2_mole_frac, _all_mole_frac, W_eq) in results_by_case.items():
        ax3.plot(W, co2_mole_frac, label=case_label)
        ax3.axvline(W_eq, color="black", linestyle="--", alpha=0.5)

    ax3.set_title("CO2 mole fraction vs reactor length")
    ax3.set_xlabel("Normalized catalyst weight W")
    ax3.set_ylabel("CO2 mole fraction [%]")
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc="best")
    fig3.tight_layout()
    fig3.savefig(os.path.join(OUTPUT_DIR, "co2_mole_percent.png"), dpi=200)

    fig4, axes4 = plt.subplots(nrows=2, ncols=1, figsize=(12, 12), sharex=True)
    for ax, (case_label, (W, _product_mole_frac, _mapped_mole_frac, _co2, all_mole_frac, W_eq)) in zip(axes4, results_by_case.items()):
        for comp, values in all_mole_frac.items():
            ax.plot(W, values, label=comp)
        ax.set_title(f"All species mole % vs reactor length ({case_label})")
        ax.set_ylabel("Mole fraction [%]")
        ax.grid(True, alpha=0.3)
        ax.axvline(W_eq, color="black", linestyle="--", alpha=0.5, label="Equilibrium")
        ax.legend(loc="best", fontsize=8, ncol=2)

    axes4[-1].set_xlabel("Normalized catalyst weight W")
    fig4.tight_layout()
    fig4.savefig(os.path.join(OUTPUT_DIR, "all_species_mole_percent.png"), dpi=200)

    try:
        plt.show()
    except Exception as exc:
        print(f"[WARN] Plot display failed: {exc}")
        print(f"[INFO] Plots saved to {OUTPUT_DIR}")


def print_final_mole_fractions(
    results_by_case: Dict[str, Tuple[np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], np.ndarray, Dict[str, np.ndarray], float]]
):
    for case_label, (W, _product_mole_frac, _mapped_mole_frac, _co2, all_mole_frac, _W_eq) in results_by_case.items():
        print("\n" + "=" * 70)
        print(f"FINAL MOLE FRACTIONS ({case_label})")
        print("=" * 70)
        for comp, values in all_mole_frac.items():
            print(f"  {comp:12s}: {values[-1]:.6f} %")
        print("=" * 70)


def main():
    cases = {
        '5% CO2 feed': 0.05,
        '10% CO2 feed': 0.10,
    }

    results_by_case: Dict[str, Tuple[np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], np.ndarray, Dict[str, np.ndarray], float]] = {}

    for label, co2_fraction in cases.items():
        inlet_flow = build_inlet_flow(co2_fraction, total_flow=1.0)
        profiles = safe_forward(SIM_CONFIG, inlet_flow)
        if profiles is None:
            raise RuntimeError(f"Solve failed for {label}")
        results_by_case[label] = profiles

    print_final_mole_fractions(results_by_case)
    plot_results(results_by_case)


if __name__ == "__main__":
    main()
