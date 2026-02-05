import logging
import os
import sys
from typing import Dict, Tuple

import numpy as np
import matplotlib.pyplot as plt
from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, value
from idaes.core import FlowsheetBlock

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "src"))

from ft_model.ft_rwgs_zeolite_reactor import FTRWGSReactor, discretize_reactor  # type: ignore[reportMissingImports]

logging.getLogger("pyomo.repn.plugins.nl_writer").setLevel(logging.ERROR)


SIM_CONFIG: Dict[str, float] = {
    # Model toggles
    'include_zeolite_reactions': True,
    'energy_balance': True,
    'pressure_drop': True,
    'ergun_pressure_drop': True,
    'heat_transfer': True,
    'mass_transfer': True,

    # Operating conditions
    'temperature': 523.15,   # K
    'pressure_bar': 20.0,    # bar
    'W_total': 1.0,          # kg catalyst
    'nfe': 10,               # spatial elements

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

    # Solver
    'max_iter': 500,
    'tol': 1e-6,
    'staged_solve': True,
}

PRODUCT_COMPONENTS = [
    'C1',
    'C2_C4',
    'C5_C12',
    'C13_plus',
    'iso_C5_C12',
    'aromatics',
    'coke',
]


def build_inlet_flow(co2_fraction: float, total_flow: float = 1.0) -> Dict[str, float]:
    co2_flow = total_flow * co2_fraction
    h2_flow = total_flow - co2_flow
    return {
        'CO2': co2_flow,
        'H2': h2_flow,
        'CO': 0.0,
        'H2O': 0.0,
        'C1': 0.0,
        'C2_C4': 0.0,
        'C5_C12': 0.0,
        'C13_plus': 0.0,
        'iso_C5_C12': 0.0,
        'aromatics': 0.0,
        'coke': 0.0,
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

    results = solver.solve(m, tee=False)
    if results.solver.termination_condition != TerminationCondition.optimal:
        return None

    return m


def extract_profiles(m) -> Tuple[np.ndarray, Dict[str, np.ndarray], np.ndarray]:
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

    co2_flow = np.array([value(reactor.flow_mol_comp[t, w, 'CO2']) for w in W_points])
    co2_mole_frac = 100.0 * co2_flow / total_flow

    return W, product_mole_frac, co2_mole_frac


def plot_results(results_by_case: Dict[str, Tuple[np.ndarray, Dict[str, np.ndarray], np.ndarray]]):
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 10), sharex=True)

    for ax, (case_label, (W, product_mole_frac, _co2)) in zip(axes, results_by_case.items()):
        for comp in PRODUCT_COMPONENTS:
            ax.plot(W, product_mole_frac[comp], label=comp)
        ax.set_title(f"Product distribution (mole %) vs reactor length ({case_label})")
        ax.set_ylabel("Mole fraction [%]")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=9)

    axes[-1].set_xlabel("Normalized catalyst weight W")
    fig.tight_layout()

    fig2, ax2 = plt.subplots(figsize=(10, 5))
    for case_label, (W, _product_mole_frac, co2_mole_frac) in results_by_case.items():
        ax2.plot(W, co2_mole_frac, label=case_label)

    ax2.set_title("CO2 mole fraction vs reactor length")
    ax2.set_xlabel("Normalized catalyst weight W")
    ax2.set_ylabel("CO2 mole fraction [%]")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best")
    fig2.tight_layout()

    plt.show()


def main():
    cases = {
        '5% CO2 feed': 0.05,
        '10% CO2 feed': 0.10,
    }

    results_by_case: Dict[str, Tuple[np.ndarray, Dict[str, np.ndarray], np.ndarray]] = {}

    for label, co2_fraction in cases.items():
        inlet_flow = build_inlet_flow(co2_fraction, total_flow=1.0)
        model = solve_case(SIM_CONFIG, inlet_flow)
        if model is None:
            raise RuntimeError(f"Solve failed for {label}")
        results_by_case[label] = extract_profiles(model)

    plot_results(results_by_case)


if __name__ == "__main__":
    main()
