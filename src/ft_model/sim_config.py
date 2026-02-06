"""
Global simulation configuration for FT/RWGS reactor runs.

This is the single source of truth for default settings used across scripts.
"""

from typing import Dict

SIM_CONFIG: Dict[str, float] = {
    # Model toggles
    "include_zeolite_reactions": True,
    "energy_balance": False,
    "pressure_drop": False,
    "ergun_pressure_drop": True,
    "heat_transfer": False,
    "mass_transfer": False,
    "kinetics_model": "marvast_2005",

    # Operating conditions
    "temperature": 523.15,   # K
    "pressure_bar": 20.0,    # bar
    "W_total": 5.0,          # kg catalyst
    "nfe": 10,               # spatial elements

    # Kinetics (base case)
    "k_rwgs": 0.001,
    "Keq_rwgs": 0.8,
    "k_c1": 0.0002,
    "k_c2_c4": 0.0001,
    "k_c5_c12": 0.00005,
    "k_c13_plus": 0.00002,
    "k_cracking": 0.0002,
    "k_light_cracking": 0.0001,
    "k_isomerization": 0.0001,
    "k_oligomerization": 0.00008,
    "k_aromatization": 0.00005,
    "k_coke_formation": 0.00001,
    "kfts_ref": 6.4e-4,
    "E_app": 23000.0,
    "b_ref": 1.6e-2,
    "dH_b": -28500.0,
    "T_ref": 543.0,

    # Transport and geometry
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

    # Selectivity/cracking factors
    "beta_gasoline": 0.7,
    "beta_jet": 0.8,
    "beta_diesel": 0.6,
    "split_c5_gasoline": 0.5,
    "split_c5_jet": 0.5,
    "split_c13_diesel": 0.7,

    # Solver
    "max_iter": 500,
    "tol": 1e-6,
    "acceptable_tol": 1e-5,
    "linear_solver": "mumps",
    "bound_push": 1e-8,
    "mu_strategy": "adaptive",
    "staged_solve": True,

    # Equilibrium detection
    "equilibrium_tol": 5e-4,
}
