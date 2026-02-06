"""
RWGS + Fischer-Tropsch + Zeolite Upgraded Packed-Bed Reactor

This reactor combines:
  1. Reverse Water-Gas Shift (RWGS) synthesis
  2. Fischer-Tropsch (FT) synthesis
  3. Lumped Zeolite Upgrading (optional)

with strict atom conservation for all reactions.

CORRECTED Reactions (atom-balanced):

FT Synthesis:
1. RWGS:    CO2 + H2   <-> CO   + H2O
2. C1:      CO  + 3H2  -> C1    + H2O
3. C2-C4:   4CO + 8H2  -> C2-C4 + 4H2O
4. C5-C12:  8CO + 16H2 -> C5-C12 + 8H2O
5. C13+:   16CO + 32H2 -> C13+  + 16H2O

Zeolite Upgrading (optional, enable with config flag):
6. Cracking:          C13+ -> 2*C5-C12
7. Light Cracking:    C5-C12 -> 2*C2-C4
8. Isomerization:     C5-C12 -> iso-C5-C12
9. Oligomerization:   2*C2-C4 -> C5-C12
10. Aromatization:    2*C2-C4 -> aromatics + 3H2
11. Coke Formation:   aromatics -> coke + H2
"""

# type: ignore  # Pyomo/IDAES type hints not fully recognized by Pylance

from typing import Dict, Optional, TYPE_CHECKING

import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, TransformationFactory, value
from pyomo.common.config import ConfigValue

from idaes.core import (
    FlowsheetBlock,
    UnitModelBlock,
    UnitModelBlockData,
    declare_process_block_class,
)
from idaes.core.util import scaling as iscale

if TYPE_CHECKING:
    class FTRWGSReactor(UnitModelBlock):
        pass


# ==================== ATOMIC COMPOSITION ====================

ATOMIC_COMPOSITION = {
    # RWGS and FT reactants/products
    'CO2': {'C': 1, 'H': 0, 'O': 2},
    'H2': {'C': 0, 'H': 2, 'O': 0},
    'CO': {'C': 1, 'H': 0, 'O': 1},
    'H2O': {'C': 0, 'H': 2, 'O': 1},
    # Hydrocarbon lumps
    'C1': {'C': 1, 'H': 4, 'O': 0},          # CH4
    'C2_C4': {'C': 4, 'H': 8, 'O': 0},           # C4H8 (olefin lump)
    'C5_C12': {'C': 6.05, 'H': 12.36, 'O': 0},   # C6.05H12.36 (gasoline-range lump)
    'C13_plus': {'C': 12.10, 'H': 24.72, 'O': 0},  # 2*C5-C12 (consistent wax lump)
    'iso_C5_C12': {'C': 6.05, 'H': 12.36, 'O': 0},
    'aromatics': {'C': 8, 'H': 10, 'O': 0},
    'coke': {'C': 1, 'H': 0, 'O': 0},
}


def check_atom_balance(reaction_name: str, stoichiometry: Dict[str, float], tolerance: float = 1e-6):
    """
    Verify atom balance for a reaction stoichiometry.
    """
    for atom in ['C', 'H', 'O']:
        total = 0.0
        for comp, coeff in stoichiometry.items():
            if comp not in ATOMIC_COMPOSITION:
                raise KeyError(f"Component '{comp}' missing atomic composition")
            total += coeff * ATOMIC_COMPOSITION[comp].get(atom, 0)

        if abs(total) > tolerance:
            raise ValueError(
                f"Reaction '{reaction_name}' has UNBALANCED {atom} atoms: {total}"
            )


def count_atoms_in_stream(flow: Dict[str, float]) -> Dict[str, float]:
    """
    Count total C/H/O atoms in a flow stream.
    """
    atoms = {'C': 0.0, 'H': 0.0, 'O': 0.0}
    for comp, mol in flow.items():
        if comp not in ATOMIC_COMPOSITION:
            raise KeyError(f"Component '{comp}' missing atomic composition")
        for atom in atoms:
            atoms[atom] += mol * ATOMIC_COMPOSITION[comp].get(atom, 0)
    return atoms


# ==================== FT REACTION STOICHIOMETRIES ====================

# RWGS: CO2 + H2 <-> CO + H2O
RWGS_STOICHIOMETRY = {
    'CO2': -1.0,
    'H2': -1.0,
    'CO': 1.0,
    'H2O': 1.0,
    'C1': 0.0,
    'C2_C4': 0.0,
    'C5_C12': 0.0,
    'C13_plus': 0.0,
    'iso_C5_C12': 0.0,
    'aromatics': 0.0,
    'coke': 0.0,
}

# C1 formation: CO + 3H2 -> CH4 + H2O
C1_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -3.0,
    'CO': -1.0,
    'H2O': 1.0,
    'C1': 1.0,
    'C2_C4': 0.0,
    'C5_C12': 0.0,
    'C13_plus': 0.0,
    'iso_C5_C12': 0.0,
    'aromatics': 0.0,
    'coke': 0.0,
}

# C2-C4 formation: 4CO + 8H2 -> C4H8 + 4H2O
C2_C4_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -8.0,
    'CO': -4.0,
    'H2O': 4.0,
    'C1': 0.0,
    'C2_C4': 1.0,
    'C5_C12': 0.0,
    'C13_plus': 0.0,
    'iso_C5_C12': 0.0,
    'aromatics': 0.0,
    'coke': 0.0,
}

# C5-C12 (gasoline-range) formation: 6.05CO + 12.23H2 -> C6.05H12.36 + 6.05H2O
C5_C12_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -12.23,
    'CO': -6.05,
    'H2O': 6.05,
    'C1': 0.0,
    'C2_C4': 0.0,
    'C5_C12': 1.0,
    'C13_plus': 0.0,
    'iso_C5_C12': 0.0,
    'aromatics': 0.0,
    'coke': 0.0,
}

# C13+ formation: 12.10CO + 24.46H2 -> C12.10H24.72 + 12.10H2O
C13_PLUS_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -24.46,
    'CO': -12.10,
    'H2O': 12.10,
    'C1': 0.0,
    'C2_C4': 0.0,
    'C5_C12': 0.0,
    'C13_plus': 1.0,
    'iso_C5_C12': 0.0,
    'aromatics': 0.0,
    'coke': 0.0,
}


# ==================== ZEOLITE REACTION STOICHIOMETRIES ====================

# Wax cracking: C12.10H24.72 -> 2*C6.05H12.36
CRACKING_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': 0.0,
    'CO': 0.0,
    'H2O': 0.0,
    'C1': 0.0,
    'C2_C4': 0.0,
    'C5_C12': 2.0,
    'C13_plus': -1.0,
    'iso_C5_C12': 0.0,
    'aromatics': 0.0,
    'coke': 0.0,
}

# Light cracking: C6.05H12.36 -> 1.5125*C4H8 + 0.13*H2
LIGHT_CRACKING_STOICHIOMETRY = {
    'CO2': 0.0,
    'CO': 0.0,
    'H2O': 0.0,
    'H2': 0.13,
    'C1': 0.0,
    'C2_C4': 1.5125,
    'C5_C12': -1.0,
    'C13_plus': 0.0,
    'iso_C5_C12': 0.0,
    'aromatics': 0.0,
    'coke': 0.0,
}

# Isomerization: C6.05H12.36 -> iso-C6.05H12.36
ISOMERIZATION_STOICHIOMETRY = {
    'CO2': 0.0,
    'CO': 0.0,
    'H2O': 0.0,
    'H2': 0.0,
    'C1': 0.0,
    'C2_C4': 0.0,
    'C5_C12': -1.0,
    'C13_plus': 0.0,
    'iso_C5_C12': 1.0,
    'aromatics': 0.0,
    'coke': 0.0,
}

# Oligomerization: 1.5125*C4H8 + 0.13*H2 -> C6.05H12.36
OLIGOMERIZATION_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -0.13,
    'CO': 0.0,
    'H2O': 0.0,
    'C1': 0.0,
    'C2_C4': -1.5125,
    'C5_C12': 1.0,
    'C13_plus': 0.0,
    'iso_C5_C12': 0.0,
    'aromatics': 0.0,
    'coke': 0.0,
}

# Aromatization: 2*C4H8 -> C8H10 + 3*H2
AROMATIZATION_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': 3.0,
    'CO': 0.0,
    'H2O': 0.0,
    'C1': 0.0,
    'C2_C4': -2.0,
    'C5_C12': 0.0,
    'C13_plus': 0.0,
    'iso_C5_C12': 0.0,
    'aromatics': 1.0,
    'coke': 0.0,
}

# Coke Formation: C8H10 -> 8*C + 5*H2
COKE_FORMATION_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': 5.0,
    'CO': 0.0,
    'H2O': 0.0,
    'C1': 0.0,
    'C2_C4': 0.0,
    'C5_C12': 0.0,
    'C13_plus': 0.0,
    'iso_C5_C12': 0.0,
    'aromatics': -1.0,
    'coke': 8.0,
}


# Dictionary of all FT reactions
FT_REACTIONS = {
    'rwgs': RWGS_STOICHIOMETRY,
    'c1': C1_STOICHIOMETRY,
    'c2_c4': C2_C4_STOICHIOMETRY,
    'c5_c12': C5_C12_STOICHIOMETRY,
    'c13_plus': C13_PLUS_STOICHIOMETRY,
}

# Dictionary of all zeolite reactions
ZEOLITE_REACTIONS = {
    'cracking': CRACKING_STOICHIOMETRY,
    'light_cracking': LIGHT_CRACKING_STOICHIOMETRY,
    'isomerization': ISOMERIZATION_STOICHIOMETRY,
    'oligomerization': OLIGOMERIZATION_STOICHIOMETRY,
    'aromatization': AROMATIZATION_STOICHIOMETRY,
    'coke_formation': COKE_FORMATION_STOICHIOMETRY,
}

# All reactions combined
ALL_REACTIONS = {
    **FT_REACTIONS,
    **ZEOLITE_REACTIONS,
}

# Verify all stoichiometries on module load
print("\n" + "="*70)
print("STOICHIOMETRY VERIFICATION")
print("="*70)
try:
    # Verify FT reactions
    check_atom_balance('RWGS', RWGS_STOICHIOMETRY)
    check_atom_balance('C1 Formation', C1_STOICHIOMETRY)
    check_atom_balance('C2-C4 Formation', C2_C4_STOICHIOMETRY)
    check_atom_balance('C5-C12 Formation', C5_C12_STOICHIOMETRY)
    check_atom_balance('C13+ Formation', C13_PLUS_STOICHIOMETRY)

    # Verify zeolite reactions
    check_atom_balance('Cracking', CRACKING_STOICHIOMETRY)
    check_atom_balance('Light Cracking', LIGHT_CRACKING_STOICHIOMETRY)
    check_atom_balance('Isomerization', ISOMERIZATION_STOICHIOMETRY)
    check_atom_balance('Oligomerization', OLIGOMERIZATION_STOICHIOMETRY)
    check_atom_balance('Aromatization', AROMATIZATION_STOICHIOMETRY)
    check_atom_balance('Coke Formation', COKE_FORMATION_STOICHIOMETRY)

    print("[OK] FT reactions verified: C/H/O atom balance OK")
    print("[OK] Zeolite reactions verified: C/H atom balance OK")
except ValueError as e:
    print(f"[ERROR] {e}")
    raise
print("="*70 + "\n")


@declare_process_block_class('FTRWGSReactor')
class FTRWGSReactorData(UnitModelBlockData):
    """
    RWGS + Fischer-Tropsch + Zeolite Packed-Bed Reactor with atom conservation.

    FT Reactions (all atom-balanced):
    1. RWGS:    CO2 + H2   <-> CO   + H2O
    2. C1:      CO  + 3H2  -> C1    + H2O
    3. C2-C4:   4CO + 8H2  -> C2-C4 + 4H2O
    4. C5-C12:  8CO + 16H2 -> C5-C12 + 8H2O
    5. C13+:   16CO + 32H2 -> C13+  + 16H2O

    Zeolite Upgrading Reactions (optional, enable via config):
    6. Cracking:          C13+ -> 2*C5-C12
    7. Light Cracking:    C5-C12 -> 2*C2-C4
    8. Isomerization:     C5-C12 -> iso-C5-C12
    9. Oligomerization:   2*C2-C4 -> C5-C12
    10. Aromatization:    2*C2-C4 -> aromatics + 3H2
    11. Coke Formation:   aromatics -> coke + H2
    """

    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        'include_zeolite_reactions',
        ConfigValue(
            default=True,
            domain=bool,
            doc='Enable zeolite upgrading reactions (default=True)',
        ),
    )
    CONFIG.declare(
        'energy_balance',
        ConfigValue(
            default=False,
            domain=bool,
            doc='Enable non-isothermal energy balance (default=False)',
        ),
    )
    CONFIG.declare(
        'pressure_drop',
        ConfigValue(
            default=False,
            domain=bool,
            doc='Enable pressure drop along catalyst bed (default=False)',
        ),
    )
    CONFIG.declare(
        'ergun_pressure_drop',
        ConfigValue(
            default=False,
            domain=bool,
            doc='Use Ergun equation for pressure drop (default=False)',
        ),
    )
    CONFIG.declare(
        'heat_transfer',
        ConfigValue(
            default=False,
            domain=bool,
            doc='Enable heat transfer to coolant/wall (default=False)',
        ),
    )
    CONFIG.declare(
        'mass_transfer',
        ConfigValue(
            default=False,
            domain=bool,
            doc='Enable effectiveness factors for mass transfer limits (default=False)',
        ),
    )
    CONFIG.declare(
        'kinetics_model',
        ConfigValue(
            default='lumped_simple',
            domain=str,
            doc='FT kinetics model: lumped_simple or marvast_2005',
        ),
    )

    def build(self):
        """
        Build the RWGS + FT + optional Zeolite reactor model.

        Creates spatial domain, state variables, parameters, and constraints.
        """
        super().build()

        # Component list - RWGS + FT lumps + zeolite products
        self.component_list = [
            'CO2', 'H2', 'CO', 'H2O',
            'C1', 'C2_C4', 'C5_C12', 'C13_plus',
            'iso_C5_C12', 'aromatics', 'coke',
        ]
        
        # FT reaction list
        self.ft_reaction_list = ['rwgs', 'c1', 'c2_c4', 'c5_c12', 'c13_plus']
        
        # Zeolite reaction list (can be disabled via config)
        self.zeo_reaction_list = [
            'cracking',
            'light_cracking',
            'isomerization',
            'oligomerization',
            'aromatization',
            'coke_formation',
        ] if self.config.include_zeolite_reactions else []

        # Spatial domain: normalized catalyst weight [0, 1]
        self.W = ContinuousSet(bounds=(0, 1.0))

        # ==================== PARAMETERS ====================

        self.W_total = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Total catalyst weight [kg]',
        )

        # RWGS reaction parameters
        self.k_rwgs = pyo.Param(
            initialize=0.1,
            mutable=True,
            doc='RWGS forward rate constant',
        )

        self.Keq_rwgs = pyo.Param(
            initialize=0.8,
            mutable=True,
            doc='RWGS equilibrium constant (dimensionless)',
        )

        # C1 formation rate constant
        self.k_c1 = pyo.Param(
            initialize=0.05,
            mutable=True,
            doc='C1 formation rate constant',
        )

        # C2-C4 formation rate constant
        self.k_c2_c4 = pyo.Param(
            initialize=0.02,
            mutable=True,
            doc='C2-C4 formation rate constant',
        )

        # C5-C12 formation rate constant
        self.k_c5_c12 = pyo.Param(
            initialize=0.01,
            mutable=True,
            doc='C5-C12 formation rate constant',
        )

        # C13+ formation rate constant
        self.k_c13_plus = pyo.Param(
            initialize=0.005,
            mutable=True,
            doc='C13+ formation rate constant',
        )

        # ==================== RWGS+FTS (2017) KINETICS PARAMETERS ====================
        # r_CDH = k_FTS(T) * P_H2 / (1 + a * P_H2O / (P_CO2 * P_H2)), a = b(T) * k_p
        self.kfts_ref = pyo.Param(
            initialize=6.4e-4,
            mutable=True,
            doc='k_FTS at reference T [mol/(g_cat*h*bar)]',
        )
        self.E_app = pyo.Param(
            initialize=23000.0,
            mutable=True,
            doc='Apparent activation energy [J/mol]',
        )
        self.b_ref = pyo.Param(
            initialize=1.6e-2,
            mutable=True,
            doc='b parameter at reference T [1/bar]',
        )
        self.dH_b = pyo.Param(
            initialize=-28500.0,
            mutable=True,
            doc='Adsorption enthalpy for b [J/mol]',
        )
        self.T_ref_kin = pyo.Param(
            initialize=543.0,
            mutable=True,
            doc='Reference temperature for k_FTS and b [K]',
        )

        # ==================== MARVAST 2005 FT KINETICS PARAMETERS ====================
        # Rj = k0 * exp(-E/RT) * P_CO^m * P_H2^n
        # Defaults from Marvast et al. (Chem Eng Technol, 2005)
        self.m_c1 = pyo.Param(initialize=-1.0889, mutable=True, doc='C1 m exponent')
        self.n_c1 = pyo.Param(initialize=1.5662, mutable=True, doc='C1 n exponent')
        self.k0_c1 = pyo.Param(initialize=142583.8, mutable=True, doc='C1 pre-exponential')
        self.E_c1 = pyo.Param(initialize=83423.9, mutable=True, doc='C1 activation energy [J/mol]')

        self.m_c2h4 = pyo.Param(initialize=0.7622, mutable=True, doc='C2H4 m exponent')
        self.n_c2h4 = pyo.Param(initialize=0.0728, mutable=True, doc='C2H4 n exponent')
        self.k0_c2h4 = pyo.Param(initialize=51.556, mutable=True, doc='C2H4 pre-exponential')
        self.E_c2h4 = pyo.Param(initialize=65018.0, mutable=True, doc='C2H4 activation energy [J/mol]')

        self.m_c2h6 = pyo.Param(initialize=-0.5645, mutable=True, doc='C2H6 m exponent')
        self.n_c2h6 = pyo.Param(initialize=1.3155, mutable=True, doc='C2H6 n exponent')
        self.k0_c2h6 = pyo.Param(initialize=24.717, mutable=True, doc='C2H6 pre-exponential')
        self.E_c2h6 = pyo.Param(initialize=49782.0, mutable=True, doc='C2H6 activation energy [J/mol]')

        self.m_c3h8 = pyo.Param(initialize=0.4051, mutable=True, doc='C3H8 m exponent')
        self.n_c3h8 = pyo.Param(initialize=0.6635, mutable=True, doc='C3H8 n exponent')
        self.k0_c3h8 = pyo.Param(initialize=0.4632, mutable=True, doc='C3H8 pre-exponential')
        self.E_c3h8 = pyo.Param(initialize=34885.5, mutable=True, doc='C3H8 activation energy [J/mol]')

        self.m_nc4h10 = pyo.Param(initialize=0.4728, mutable=True, doc='n-C4H10 m exponent')
        self.n_nc4h10 = pyo.Param(initialize=1.1389, mutable=True, doc='n-C4H10 n exponent')
        self.k0_nc4h10 = pyo.Param(initialize=0.00474, mutable=True, doc='n-C4H10 pre-exponential')
        self.E_nc4h10 = pyo.Param(initialize=27728.9, mutable=True, doc='n-C4H10 activation energy [J/mol]')

        self.m_ic4h10 = pyo.Param(initialize=0.8204, mutable=True, doc='i-C4H10 m exponent')
        self.n_ic4h10 = pyo.Param(initialize=0.5026, mutable=True, doc='i-C4H10 n exponent')
        self.k0_ic4h10 = pyo.Param(initialize=0.00832, mutable=True, doc='i-C4H10 pre-exponential')
        self.E_ic4h10 = pyo.Param(initialize=25730.1, mutable=True, doc='i-C4H10 activation energy [J/mol]')

        self.m_c5_c12 = pyo.Param(initialize=0.5850, mutable=True, doc='C5-C12 m exponent (C6.05H12.36)')
        self.n_c5_c12 = pyo.Param(initialize=0.5982, mutable=True, doc='C5-C12 n exponent (C6.05H12.36)')
        self.k0_c5_c12 = pyo.Param(initialize=0.02316, mutable=True, doc='C5-C12 pre-exponential (C6.05H12.36)')
        self.E_c5_c12 = pyo.Param(initialize=23564.3, mutable=True, doc='C5-C12 activation energy [J/mol]')

        self.dp_dw = pyo.Param(
            initialize=0.0,
            mutable=True,
            doc='Pressure drop per kg catalyst [Pa/kg]',
        )

        self.ergun_porosity = pyo.Param(
            initialize=0.40,
            mutable=True,
            doc='Bed void fraction for Ergun equation [-]',
        )
        self.particle_diameter = pyo.Param(
            initialize=2.0e-3,
            mutable=True,
            doc='Particle diameter [m]',
        )
        self.catalyst_bulk_density = pyo.Param(
            initialize=1200.0,
            mutable=True,
            doc='Catalyst bulk density [kg/m^3]',
        )
        self.reactor_diameter = pyo.Param(
            initialize=0.10,
            mutable=True,
            doc='Reactor inner diameter [m]',
        )
        self.reactor_length = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Reactor length [m]',
        )
        self.gas_viscosity = pyo.Param(
            initialize=2.0e-5,
            mutable=True,
            doc='Gas viscosity [Pa*s]',
        )

        self.ua_per_kg = pyo.Param(
            initialize=0.0,
            mutable=True,
            doc='Overall heat transfer coefficient per kg catalyst [kJ/(s*kg*K)]',
        )
        self.T_coolant = pyo.Param(
            initialize=500.0,
            mutable=True,
            doc='Coolant/wall temperature [K]',
        )

        self.eta_ft = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Effectiveness factor for FT/RWGS reactions [-]',
        )
        self.eta_zeolite = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Effectiveness factor for zeolite reactions [-]',
        )

        self.R_gas = pyo.Param(
            initialize=8.314e3,
            mutable=True,
            doc='Gas constant [Pa*m^3/(kmol*K)]',
        )

        self.mw_comp = pyo.Param(
            self.component_list,
            initialize={
                'CO2': 44.01,
                'H2': 2.016,
                'CO': 28.01,
                'H2O': 18.015,
                'C1': 16.04,
                'C2_C4': 56.11,
                'C5_C12': 112.21,
                'C13_plus': 224.43,
                'iso_C5_C12': 112.21,
                'aromatics': 106.17,
                'coke': 12.01,
            },
            mutable=True,
            doc='Component molecular weights [kg/kmol]',
        )

        # ==================== ZEOLITE REACTION PARAMETERS ====================

        self.k_cracking = pyo.Param(
            initialize=0.05,
            mutable=True,
            doc='Wax cracking rate constant [kmol/(kg_cat·s)]',
        )

        self.k_light_cracking = pyo.Param(
            initialize=0.03,
            mutable=True,
            doc='Light cracking rate constant [kmol/(kg_cat·s)]',
        )

        self.k_isomerization = pyo.Param(
            initialize=0.02,
            mutable=True,
            doc='Isomerization rate constant [kmol/(kg_cat·s)]',
        )

        self.k_oligomerization = pyo.Param(
            initialize=0.015,
            mutable=True,
            doc='Oligomerization rate constant [kmol/(kg_cat·s)]',
        )

        self.k_aromatization = pyo.Param(
            initialize=0.01,
            mutable=True,
            doc='Aromatization rate constant [kmol/(kg_cat·s)]',
        )

        self.k_coke_formation = pyo.Param(
            initialize=0.001,
            mutable=True,
            doc='Coke formation rate constant [kmol/(kg_cat·s)]',
        )

        self.beta_gasoline = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Cracking degree for gasoline-range [-]',
        )
        self.beta_jet = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Cracking degree for jet-range [-]',
        )
        self.beta_diesel = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Cracking degree for diesel-range [-]',
        )
        self.split_c5_gasoline = pyo.Param(
            initialize=0.5,
            mutable=True,
            doc='Fraction of C5_C12 treated as gasoline-range [-]',
        )
        self.split_c5_jet = pyo.Param(
            initialize=0.5,
            mutable=True,
            doc='Fraction of C5_C12 treated as jet-range [-]',
        )
        self.split_c13_diesel = pyo.Param(
            initialize=0.7,
            mutable=True,
            doc='Fraction of C13_plus treated as diesel-range [-]',
        )

        # ==================== STATE VARIABLES ====================

        self.flow_mol_comp = pyo.Var(
            self.flowsheet().time,
            self.W,
            self.component_list,
            initialize=0.1,
            bounds=(0.0, None),
            doc='Component molar flow rates [kmol/s]',
        )

        self.temperature = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=500.0,
            bounds=(200.0, 1000.0),
            doc='Temperature [K]',
        )

        self.pressure = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=2e6,
            bounds=(1e5, 5e6),
            doc='Pressure [Pa]',
        )

        self.flow_mol_total = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=1.0,
            bounds=(1e-6, None),
            doc='Total molar flow [kmol/s]',
        )

        self.mole_frac = pyo.Expression(
            self.flowsheet().time,
            self.W,
            self.component_list,
            rule=lambda b, t, w, c: b.flow_mol_comp[t, w, c] / (b.flow_mol_total[t, w] + 1e-12),
            doc='Mole fractions',
        )

        self.partial_pressure = pyo.Expression(
            self.flowsheet().time,
            self.W,
            self.component_list,
            rule=lambda b, t, w, c: b.mole_frac[t, w, c] * b.pressure[t, w],
            doc='Partial pressures [Pa]',
        )

        self.flow_gasoline = pyo.Expression(
            self.flowsheet().time,
            self.W,
            rule=lambda b, t, w: b.split_c5_gasoline * b.flow_mol_comp[t, w, 'C5_C12']
            + b.flow_mol_comp[t, w, 'iso_C5_C12'],
            doc='Gasoline-range flow [kmol/s] (mapped from C5_C12 and iso_C5_C12)',
        )
        self.flow_jet = pyo.Expression(
            self.flowsheet().time,
            self.W,
            rule=lambda b, t, w: b.split_c5_jet * b.flow_mol_comp[t, w, 'C5_C12'],
            doc='Jet-range flow [kmol/s] (mapped from C5_C12)',
        )
        self.flow_diesel = pyo.Expression(
            self.flowsheet().time,
            self.W,
            rule=lambda b, t, w: b.split_c13_diesel * b.flow_mol_comp[t, w, 'C13_plus'],
            doc='Diesel-range flow [kmol/s] (mapped from C13_plus)',
        )

        # ==================== THERMODYNAMICS (PLACEHOLDER PROPERTY BLOCK) ====================

        self.T_ref = pyo.Param(
            initialize=298.15,
            mutable=True,
            doc='Reference temperature for enthalpy [K]',
        )

        # Constant heat capacities [kJ/kmol-K] (placeholder values)
        self.cp_comp = pyo.Param(
            self.component_list,
            initialize={
                'CO2': 37.0,
                'H2': 29.0,
                'CO': 29.0,
                'H2O': 34.0,
                'C1': 35.0,
                'C2_C4': 95.0,
                'C5_C12': 180.0,
                'C13_plus': 240.0,
                'iso_C5_C12': 180.0,
                'aromatics': 140.0,
                'coke': 8.0,
            },
            mutable=True,
            doc='Component heat capacities [kJ/kmol-K]',
        )

        # Standard enthalpies of formation [kJ/kmol] (placeholder values)
        self.h_form = pyo.Param(
            self.component_list,
            initialize={
                'CO2': -393520.0,
                'H2': 0.0,
                'CO': -110530.0,
                'H2O': -241820.0,
                'C1': -74850.0,
                'C2_C4': -20000.0,
                'C5_C12': -120000.0,
                'C13_plus': -220000.0,
                'iso_C5_C12': -120000.0,
                'aromatics': 83000.0,
                'coke': 0.0,
            },
            mutable=True,
            doc='Enthalpy of formation [kJ/kmol]',
        )

        self.h_comp = pyo.Expression(
            self.flowsheet().time,
            self.W,
            self.component_list,
            doc='Component molar enthalpy [kJ/kmol]',
            rule=lambda b, t, w, c: b.h_form[c] + b.cp_comp[c] * (b.temperature[t, w] - b.T_ref),
        )

        def _props_rule(b, t, w):
            b.enth_mol = pyo.Expression(
                expr=sum(
                    self.mole_frac[t, w, c] * self.h_comp[t, w, c]
                    for c in self.component_list
                )
            )
            b.cp_mol = pyo.Expression(
                expr=sum(
                    self.mole_frac[t, w, c] * self.cp_comp[c]
                    for c in self.component_list
                )
            )

        self.props = pyo.Block(self.flowsheet().time, self.W, rule=_props_rule)

        # ==================== TRANSPORT/GEOMETRY EXPRESSIONS ====================

        self.area = pyo.Expression(
            expr=3.141592653589793 * (self.reactor_diameter / 2.0) ** 2,
            doc='Reactor cross-sectional area [m^2]',
        )

        self.mw_mix = pyo.Expression(
            self.flowsheet().time,
            self.W,
            rule=lambda b, t, w: sum(
                b.mole_frac[t, w, c] * b.mw_comp[c] for c in b.component_list
            ),
            doc='Mixture molecular weight [kg/kmol]',
        )

        self.gas_density = pyo.Expression(
            self.flowsheet().time,
            self.W,
            rule=lambda b, t, w: b.pressure[t, w] * b.mw_mix[t, w] / (b.R_gas * b.temperature[t, w]),
            doc='Ideal-gas density [kg/m^3]',
        )

        self.vol_flow = pyo.Expression(
            self.flowsheet().time,
            self.W,
            rule=lambda b, t, w: b.flow_mol_total[t, w] * b.R_gas * b.temperature[t, w] / b.pressure[t, w],
            doc='Volumetric flow [m^3/s]',
        )

        self.superficial_velocity = pyo.Expression(
            self.flowsheet().time,
            self.W,
            rule=lambda b, t, w: b.vol_flow[t, w] / b.area,
            doc='Superficial velocity [m/s]',
        )

        self.dP_dz_ergun = pyo.Expression(
            self.flowsheet().time,
            self.W,
            rule=lambda b, t, w: (
                150.0 * (1.0 - b.ergun_porosity) ** 2 * b.gas_viscosity
                * b.superficial_velocity[t, w]
                / (b.particle_diameter ** 2 * b.ergun_porosity ** 3)
                + 1.75 * (1.0 - b.ergun_porosity)
                * b.gas_density[t, w] * b.superficial_velocity[t, w] ** 2
                / (b.particle_diameter * b.ergun_porosity ** 3)
            ),
            doc='Ergun pressure gradient [Pa/m]',
        )

        self.dP_dW_ergun = pyo.Expression(
            self.flowsheet().time,
            self.W,
            rule=lambda b, t, w: (
                b.dP_dz_ergun[t, w]
                / ((1.0 - b.ergun_porosity) * b.catalyst_bulk_density * b.area)
                * b.W_total
            ),
            doc='Ergun pressure gradient in normalized W [Pa]',
        )

        self.rate_multiplier = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Homotopy multiplier for reaction rates (0 to 1)',
        )

        self.pressure_drop_multiplier = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Homotopy multiplier for pressure drop (0 to 1)',
        )

        # ==================== REACTION RATES ====================
        self.rate_rwgs = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.01,
            bounds=None,
            doc='RWGS reaction rate [kmol/(kg_cat·s)]',
        )
        
        self.rate_c1 = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.001,
            bounds=None,
            doc='C1 formation rate [kmol/(kg_cat·s)]',
        )
        
        self.rate_c2_c4 = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.0005,
            bounds=None,
            doc='C2-C4 formation rate [kmol/(kg_cat·s)]',
        )
        
        self.rate_c5_c12 = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.0001,
            bounds=None,
            doc='C5-C12 formation rate [kmol/(kg_cat·s)]',
        )

        self.rate_c13_plus = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.00005,
            bounds=None,
            doc='C13+ formation rate [kmol/(kg_cat·s)]',
        )
        
        # ==================== ZEOLITE REACTION RATES ====================
        
        self.rate_cracking = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.0001 if self.config.include_zeolite_reactions else 0.0,
            bounds=None,
            doc='Wax cracking rate [kmol/(kg_cat·s)]',
        )
        
        self.rate_light_cracking = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.00005 if self.config.include_zeolite_reactions else 0.0,
            bounds=None,
            doc='Light cracking rate [kmol/(kg_cat·s)]',
        )
        
        self.rate_isomerization = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.00001 if self.config.include_zeolite_reactions else 0.0,
            bounds=None,
            doc='Isomerization rate [kmol/(kg_cat·s)]',
        )

        self.rate_oligomerization = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.00001 if self.config.include_zeolite_reactions else 0.0,
            bounds=None,
            doc='Oligomerization rate [kmol/(kg_cat·s)]',
        )

        self.rate_aromatization = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.000005 if self.config.include_zeolite_reactions else 0.0,
            bounds=None,
            doc='Aromatization rate [kmol/(kg_cat·s)]',
        )
        
        self.rate_coke_formation = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.000001 if self.config.include_zeolite_reactions else 0.0,
            bounds=None,
            doc='Coke formation rate [kmol/(kg_cat·s)]',
        )
        
        # ==================== DERIVATIVES ====================
        
        # Material balance derivatives dF/dW
        self.dF_dW = DerivativeVar(
            self.flow_mol_comp,
            wrt=self.W,
            doc='Derivative of flow w.r.t. catalyst weight',
        )

        # Temperature derivative dT/dW
        self.dT_dW = DerivativeVar(
            self.temperature,
            wrt=self.W,
            doc='Derivative of temperature w.r.t. catalyst weight',
        )

        # Pressure derivative dP/dW (used if pressure_drop is enabled)
        self.dP_dW = DerivativeVar(
            self.pressure,
            wrt=self.W,
            doc='Derivative of pressure w.r.t. catalyst weight',
        )
        
        # ==================== CONSTRAINTS ====================
        
        # Total flow calculation
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='Total molar flow',
        )
        def total_flow_eq(b, t, w):
            return b.flow_mol_total[t, w] == sum(
                b.flow_mol_comp[t, w, c] for c in b.component_list
            )
        
        if self.config.pressure_drop:
            @self.Constraint(
                self.flowsheet().time,
                self.W,
                doc='Pressure drop along catalyst bed',
            )
            def pressure_drop_eq(b, t, w):
                if b.config.ergun_pressure_drop:
                    return b.dP_dW[t, w] == -b.pressure_drop_multiplier * b.dP_dW_ergun[t, w]
                return b.dP_dW[t, w] == -b.pressure_drop_multiplier * b.dp_dw * b.W_total
        
        # RWGS rate: r = k * (p_CO2 * p_H2 - p_CO * p_H2O / Keq)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='RWGS reaction rate',
        )
        def rate_rwgs_eq(b, t, w):
            p_CO2 = b.partial_pressure[t, w, 'CO2'] / 1e5 + 1e-6
            p_H2 = b.partial_pressure[t, w, 'H2'] / 1e5 + 1e-6
            p_CO = b.partial_pressure[t, w, 'CO'] / 1e5 + 1e-6
            p_H2O = b.partial_pressure[t, w, 'H2O'] / 1e5 + 1e-6
            
            forward = p_CO2 * p_H2
            reverse = p_CO * p_H2O / b.Keq_rwgs
            
            return b.rate_rwgs[t, w] == b.rate_multiplier * b.k_rwgs * (forward - reverse)
        
        # C1 formation rate: r = k * p_CO * p_H2
        # Simplified power-law for numerical stability
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='C1 formation rate (simplified power-law)',
        )
        def rate_c1_eq(b, t, w):
            p_CO = b.partial_pressure[t, w, 'CO'] / 1e5 + 1e-6
            p_H2 = b.partial_pressure[t, w, 'H2'] / 1e5 + 1e-6
            if b.config.kinetics_model == 'marvast_2005':
                p_CO_safe = p_CO + 1e-6
                p_H2_safe = p_H2 + 1e-6
                k_T = b.k0_c1 * pyo.exp(-b.E_c1 / (8.314 * b.temperature[t, w]))
                return b.rate_c1[t, w] == b.rate_multiplier * k_T * p_CO_safe**b.m_c1 * p_H2_safe**b.n_c1
            if b.config.kinetics_model == 'rwgs_2017':
                return b.rate_c1[t, w] == 0.0
            return b.rate_c1[t, w] == b.rate_multiplier * b.k_c1 * p_CO * p_H2
        
        # C2-C4 formation rate: r = k * p_CO * p_H2
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='C2-C4 formation rate (simplified power-law)',
        )
        def rate_c2_c4_eq(b, t, w):
            p_CO = b.partial_pressure[t, w, 'CO'] / 1e5 + 1e-6
            p_H2 = b.partial_pressure[t, w, 'H2'] / 1e5 + 1e-6
            if b.config.kinetics_model == 'marvast_2005':
                p_CO_safe = p_CO + 1e-6
                p_H2_safe = p_H2 + 1e-6
                k_c2h4 = b.k0_c2h4 * pyo.exp(-b.E_c2h4 / (8.314 * b.temperature[t, w]))
                k_c2h6 = b.k0_c2h6 * pyo.exp(-b.E_c2h6 / (8.314 * b.temperature[t, w]))
                k_c3h8 = b.k0_c3h8 * pyo.exp(-b.E_c3h8 / (8.314 * b.temperature[t, w]))
                k_nc4 = b.k0_nc4h10 * pyo.exp(-b.E_nc4h10 / (8.314 * b.temperature[t, w]))
                k_ic4 = b.k0_ic4h10 * pyo.exp(-b.E_ic4h10 / (8.314 * b.temperature[t, w]))
                rate_sum = (
                    k_c2h4 * p_CO_safe**b.m_c2h4 * p_H2_safe**b.n_c2h4
                    + k_c2h6 * p_CO_safe**b.m_c2h6 * p_H2_safe**b.n_c2h6
                    + k_c3h8 * p_CO_safe**b.m_c3h8 * p_H2_safe**b.n_c3h8
                    + k_nc4 * p_CO_safe**b.m_nc4h10 * p_H2_safe**b.n_nc4h10
                    + k_ic4 * p_CO_safe**b.m_ic4h10 * p_H2_safe**b.n_ic4h10
                )
                return b.rate_c2_c4[t, w] == b.rate_multiplier * rate_sum
            if b.config.kinetics_model == 'rwgs_2017':
                return b.rate_c2_c4[t, w] == 0.0
            return b.rate_c2_c4[t, w] == b.rate_multiplier * b.k_c2_c4 * p_CO * p_H2
        
        # C5-C12 formation rate: r = k * p_CO * p_H2
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='C5-C12 formation rate (simplified power-law)',
        )
        def rate_c5_c12_eq(b, t, w):
            p_CO = b.partial_pressure[t, w, 'CO'] / 1e5 + 1e-6
            p_H2 = b.partial_pressure[t, w, 'H2'] / 1e5 + 1e-6
            p_CO2 = b.partial_pressure[t, w, 'CO2'] / 1e5 + 1e-6
            p_H2O = b.partial_pressure[t, w, 'H2O'] / 1e5 + 1e-6
            if b.config.kinetics_model == 'marvast_2005':
                p_CO_safe = p_CO + 1e-6
                p_H2_safe = p_H2 + 1e-6
                k_T = b.k0_c5_c12 * pyo.exp(-b.E_c5_c12 / (8.314 * b.temperature[t, w]))
                return b.rate_c5_c12[t, w] == b.rate_multiplier * k_T * p_CO_safe**b.m_c5_c12 * p_H2_safe**b.n_c5_c12
            if b.config.kinetics_model == 'rwgs_2017':
                return b.rate_c5_c12[t, w] == 0.0
            return b.rate_c5_c12[t, w] == b.rate_multiplier * b.k_c5_c12 * p_CO * p_H2

        # C13+ formation rate: r = k * p_CO * p_H2
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='C13+ formation rate (simplified power-law)',
        )
        def rate_c13_plus_eq(b, t, w):
            p_CO = b.partial_pressure[t, w, 'CO'] / 1e5 + 1e-6
            p_H2 = b.partial_pressure[t, w, 'H2'] / 1e5 + 1e-6
            if b.config.kinetics_model == 'marvast_2005':
                return b.rate_c13_plus[t, w] == 0.0
            if b.config.kinetics_model == 'rwgs_2017':
                return b.rate_c13_plus[t, w] == 0.0
            return b.rate_c13_plus[t, w] == b.rate_multiplier * b.k_c13_plus * p_CO * p_H2
        
        # ==================== ZEOLITE RATE EXPRESSIONS ====================
        # All zeolite reactions use simple first-order kinetics in key reactants
        
        if self.config.include_zeolite_reactions:
            # Wax cracking: r = k * C13_plus mole fraction
            @self.Constraint(
                self.flowsheet().time,
                self.W,
                doc='Wax cracking rate (first-order in C13_plus)',
            )
            def rate_cracking_eq(b, t, w):
                y_c13 = b.mole_frac[t, w, 'C13_plus']
                return b.rate_cracking[t, w] == b.rate_multiplier * b.k_cracking * b.beta_diesel * y_c13

            # Light cracking: r = k * C5_C12 mole fraction
            @self.Constraint(
                self.flowsheet().time,
                self.W,
                doc='Light cracking rate (first-order in C5_C12)',
            )
            def rate_light_cracking_eq(b, t, w):
                y_c5 = b.mole_frac[t, w, 'C5_C12']
                beta_c5 = b.beta_gasoline * b.split_c5_gasoline + b.beta_jet * b.split_c5_jet
                return b.rate_light_cracking[t, w] == b.rate_multiplier * b.k_light_cracking * beta_c5 * y_c5

            # Isomerization: r = k * C5_C12 mole fraction
            @self.Constraint(
                self.flowsheet().time,
                self.W,
                doc='Isomerization rate (first-order in C5_C12)',
            )
            def rate_isomerization_eq(b, t, w):
                y_c5 = b.mole_frac[t, w, 'C5_C12']
                return b.rate_isomerization[t, w] == b.rate_multiplier * b.k_isomerization * y_c5

            # Oligomerization: r = k * C2_C4 mole fraction
            @self.Constraint(
                self.flowsheet().time,
                self.W,
                doc='Oligomerization rate (first-order in C2_C4)',
            )
            def rate_oligomerization_eq(b, t, w):
                y_c2_c4 = b.mole_frac[t, w, 'C2_C4']
                return b.rate_oligomerization[t, w] == b.rate_multiplier * b.k_oligomerization * y_c2_c4

            # Aromatization: r = k * C2_C4 mole fraction
            @self.Constraint(
                self.flowsheet().time,
                self.W,
                doc='Aromatization rate (first-order in C2_C4)',
            )
            def rate_aromatization_eq(b, t, w):
                y_c2_c4 = b.mole_frac[t, w, 'C2_C4']
                return b.rate_aromatization[t, w] == b.rate_multiplier * b.k_aromatization * y_c2_c4

            # Coke formation: r = k * aromatics partial pressure
            @self.Constraint(
                self.flowsheet().time,
                self.W,
                doc='Coke formation rate (first-order in aromatics)',
            )
            def rate_coke_formation_eq(b, t, w):
                y_arom = b.mole_frac[t, w, 'aromatics']
                return b.rate_coke_formation[t, w] == b.rate_multiplier * b.k_coke_formation * y_arom
        else:
            @self.Constraint(self.flowsheet().time, self.W, doc='Cracking disabled')
            def rate_cracking_off(b, t, w):
                return b.rate_cracking[t, w] == 0.0

            @self.Constraint(self.flowsheet().time, self.W, doc='Light cracking disabled')
            def rate_light_cracking_off(b, t, w):
                return b.rate_light_cracking[t, w] == 0.0

            @self.Constraint(self.flowsheet().time, self.W, doc='Isomerization disabled')
            def rate_isomerization_off(b, t, w):
                return b.rate_isomerization[t, w] == 0.0

            @self.Constraint(self.flowsheet().time, self.W, doc='Oligomerization disabled')
            def rate_oligomerization_off(b, t, w):
                return b.rate_oligomerization[t, w] == 0.0

            @self.Constraint(self.flowsheet().time, self.W, doc='Aromatization disabled')
            def rate_aromatization_off(b, t, w):
                return b.rate_aromatization[t, w] == 0.0

            @self.Constraint(self.flowsheet().time, self.W, doc='Coke formation disabled')
            def rate_coke_formation_off(b, t, w):
                return b.rate_coke_formation[t, w] == 0.0

        # ==================== MATERIAL BALANCES ====================

        @self.Constraint(
            self.flowsheet().time,
            self.W,
            self.component_list,
            doc='Component material balances',
        )
        def material_balance_eq(b, t, w, c):
            ft_term = (
                RWGS_STOICHIOMETRY.get(c, 0.0) * b.rate_rwgs[t, w]
                + C1_STOICHIOMETRY.get(c, 0.0) * b.rate_c1[t, w]
                + C2_C4_STOICHIOMETRY.get(c, 0.0) * b.rate_c2_c4[t, w]
                + C5_C12_STOICHIOMETRY.get(c, 0.0) * b.rate_c5_c12[t, w]
                + C13_PLUS_STOICHIOMETRY.get(c, 0.0) * b.rate_c13_plus[t, w]
            )

            zeo_term = (
                CRACKING_STOICHIOMETRY.get(c, 0.0) * b.rate_cracking[t, w]
                + LIGHT_CRACKING_STOICHIOMETRY.get(c, 0.0) * b.rate_light_cracking[t, w]
                + ISOMERIZATION_STOICHIOMETRY.get(c, 0.0) * b.rate_isomerization[t, w]
                + OLIGOMERIZATION_STOICHIOMETRY.get(c, 0.0) * b.rate_oligomerization[t, w]
                + AROMATIZATION_STOICHIOMETRY.get(c, 0.0) * b.rate_aromatization[t, w]
                + COKE_FORMATION_STOICHIOMETRY.get(c, 0.0) * b.rate_coke_formation[t, w]
            )
            eta_ft = b.eta_ft if b.config.mass_transfer else 1.0
            eta_zeo = b.eta_zeolite if b.config.mass_transfer else 1.0

            return b.dF_dW[t, w, c] == b.W_total * (eta_ft * ft_term + eta_zeo * zeo_term)

        # ==================== REACTION ENTHALPIES ====================

        self.reaction_stoich = {
            'rwgs': RWGS_STOICHIOMETRY,
            'c1': C1_STOICHIOMETRY,
            'c2_c4': C2_C4_STOICHIOMETRY,
            'c5_c12': C5_C12_STOICHIOMETRY,
            'c13_plus': C13_PLUS_STOICHIOMETRY,
            'cracking': CRACKING_STOICHIOMETRY,
            'light_cracking': LIGHT_CRACKING_STOICHIOMETRY,
            'isomerization': ISOMERIZATION_STOICHIOMETRY,
            'oligomerization': OLIGOMERIZATION_STOICHIOMETRY,
            'aromatization': AROMATIZATION_STOICHIOMETRY,
            'coke_formation': COKE_FORMATION_STOICHIOMETRY,
        }

        self.reaction_rate = {
            'rwgs': self.rate_rwgs,
            'c1': self.rate_c1,
            'c2_c4': self.rate_c2_c4,
            'c5_c12': self.rate_c5_c12,
            'c13_plus': self.rate_c13_plus,
            'cracking': self.rate_cracking,
            'light_cracking': self.rate_light_cracking,
            'isomerization': self.rate_isomerization,
            'oligomerization': self.rate_oligomerization,
            'aromatization': self.rate_aromatization,
            'coke_formation': self.rate_coke_formation,
        }

        self.dH_rxn = pyo.Expression(
            self.flowsheet().time,
            self.W,
            self.reaction_stoich.keys(),
            doc='Reaction enthalpies [kJ/kmol]',
            rule=lambda b, t, w, r: sum(
                self.reaction_stoich[r].get(c, 0.0) * b.h_comp[t, w, c]
                for c in self.component_list
            ),
        )

        # ==================== ENERGY BALANCE ====================

        if self.config.energy_balance:
            @self.Constraint(
                self.flowsheet().time,
                self.W,
                doc='Energy balance (adiabatic)',
            )
            def energy_balance_eq(b, t, w):
                lhs = sum(
                    b.dF_dW[t, w, c] * b.h_comp[t, w, c]
                    for c in b.component_list
                ) + b.flow_mol_total[t, w] * b.props[t, w].cp_mol * b.dT_dW[t, w]

                rhs = -b.W_total * sum(
                    b.dH_rxn[t, w, r] * b.reaction_rate[r][t, w]
                    for r in b.reaction_stoich.keys()
                )

                if b.config.heat_transfer:
                    rhs += -b.ua_per_kg * b.W_total * (b.temperature[t, w] - b.T_coolant)

                return lhs == rhs
        
    def initialize(
        self,
        inlet_flow: Dict[str, float],
        temperature: float = 523.15,
        pressure: float = 2e6,
        W_total: float = 1.0,
        k_rwgs: Optional[float] = None,
        Keq_rwgs: Optional[float] = None,
        k_c1: Optional[float] = None,
        k_c2_c4: Optional[float] = None,
        k_c5_c12: Optional[float] = None,
        k_c13_plus: Optional[float] = None,
        k_cracking: Optional[float] = None,
        k_light_cracking: Optional[float] = None,
        k_isomerization: Optional[float] = None,
        k_oligomerization: Optional[float] = None,
        k_aromatization: Optional[float] = None,
        k_coke_formation: Optional[float] = None,
        kfts_ref: Optional[float] = None,
        E_app: Optional[float] = None,
        b_ref: Optional[float] = None,
        dH_b: Optional[float] = None,
        T_ref: Optional[float] = None,
        beta_gasoline: Optional[float] = None,
        beta_jet: Optional[float] = None,
        beta_diesel: Optional[float] = None,
        split_c5_gasoline: Optional[float] = None,
        split_c5_jet: Optional[float] = None,
        split_c13_diesel: Optional[float] = None,
        dp_dw: Optional[float] = None,
        ergun_porosity: Optional[float] = None,
        particle_diameter: Optional[float] = None,
        catalyst_bulk_density: Optional[float] = None,
        reactor_diameter: Optional[float] = None,
        reactor_length: Optional[float] = None,
        gas_viscosity: Optional[float] = None,
        ua_per_kg: Optional[float] = None,
        T_coolant: Optional[float] = None,
        eta_ft: Optional[float] = None,
        eta_zeolite: Optional[float] = None,
    ):
        """
        Initialize the reactor with inlet conditions and optional parameters.
        """
        t = self.flowsheet().time.first()
        w0 = self.W.first()

        # Update parameters if provided
        if W_total is not None:
            self.W_total.set_value(W_total)
        if k_rwgs is not None:
            self.k_rwgs.set_value(k_rwgs)
        if Keq_rwgs is not None:
            self.Keq_rwgs.set_value(Keq_rwgs)
        if k_c1 is not None:
            self.k_c1.set_value(k_c1)
        if k_c2_c4 is not None:
            self.k_c2_c4.set_value(k_c2_c4)
        if k_c5_c12 is not None:
            self.k_c5_c12.set_value(k_c5_c12)
        if k_c13_plus is not None:
            self.k_c13_plus.set_value(k_c13_plus)
        if k_cracking is not None:
            self.k_cracking.set_value(k_cracking)
        if k_light_cracking is not None:
            self.k_light_cracking.set_value(k_light_cracking)
        if k_isomerization is not None:
            self.k_isomerization.set_value(k_isomerization)
        if k_oligomerization is not None:
            self.k_oligomerization.set_value(k_oligomerization)
        if k_aromatization is not None:
            self.k_aromatization.set_value(k_aromatization)
        if k_coke_formation is not None:
            self.k_coke_formation.set_value(k_coke_formation)
        if kfts_ref is not None:
            self.kfts_ref.set_value(kfts_ref)
        if E_app is not None:
            self.E_app.set_value(E_app)
        if b_ref is not None:
            self.b_ref.set_value(b_ref)
        if dH_b is not None:
            self.dH_b.set_value(dH_b)
        if T_ref is not None:
            self.T_ref_kin.set_value(T_ref)
        if beta_gasoline is not None:
            self.beta_gasoline.set_value(beta_gasoline)
        if beta_jet is not None:
            self.beta_jet.set_value(beta_jet)
        if beta_diesel is not None:
            self.beta_diesel.set_value(beta_diesel)
        if split_c5_gasoline is not None:
            self.split_c5_gasoline.set_value(split_c5_gasoline)
        if split_c5_jet is not None:
            self.split_c5_jet.set_value(split_c5_jet)
        if split_c13_diesel is not None:
            self.split_c13_diesel.set_value(split_c13_diesel)
        if dp_dw is not None:
            self.dp_dw.set_value(dp_dw)
        if ergun_porosity is not None:
            self.ergun_porosity.set_value(ergun_porosity)
        if particle_diameter is not None:
            self.particle_diameter.set_value(particle_diameter)
        if catalyst_bulk_density is not None:
            self.catalyst_bulk_density.set_value(catalyst_bulk_density)
        if reactor_diameter is not None:
            self.reactor_diameter.set_value(reactor_diameter)
        if reactor_length is not None:
            self.reactor_length.set_value(reactor_length)
        if gas_viscosity is not None:
            self.gas_viscosity.set_value(gas_viscosity)
        if ua_per_kg is not None:
            self.ua_per_kg.set_value(ua_per_kg)
        if T_coolant is not None:
            self.T_coolant.set_value(T_coolant)
        if eta_ft is not None:
            self.eta_ft.set_value(eta_ft)
        if eta_zeolite is not None:
            self.eta_zeolite.set_value(eta_zeolite)

        # Fix inlet conditions at W=0
        for comp in self.component_list:
            flow = inlet_flow.get(comp, 0.0)
            self.flow_mol_comp[t, w0, comp].fix(flow)

        if self.config.pressure_drop:
            # Pressure drop enabled: fix inlet only and initialize profile
            self.pressure[t, w0].fix(pressure)
            for w in self.W:
                if w != w0:
                    self.pressure[t, w].set_value(
                        pressure - self.dp_dw.value * self.W_total.value * w
                    )
                    self.pressure[t, w].unfix()
        else:
            # Isobaric: fix pressure along W
            for w in self.W:
                self.pressure[t, w].fix(pressure)

        if self.config.energy_balance:
            # Energy balance enabled: fix inlet T only, initialize profile
            self.temperature[t, w0].fix(temperature)
            for w in self.W:
                if w != w0:
                    self.temperature[t, w].set_value(temperature)
                    self.temperature[t, w].unfix()
        else:
            # Isothermal mode: fix T along W
            for w in self.W:
                self.temperature[t, w].fix(temperature)

        # Propagate initial guesses along W
        for w in self.W:
            if w != w0:
                for comp in self.component_list:
                    self.flow_mol_comp[t, w, comp].set_value(
                        self.flow_mol_comp[t, w0, comp].value
                    )

        for w in self.W:
            total_flow = sum(
                self.flow_mol_comp[t, w, comp].value
                for comp in self.component_list
            )
            self.flow_mol_total[t, w].set_value(total_flow)

        rate_mult = self.rate_multiplier.value
        eta_ft = self.eta_ft.value
        eta_zeo = self.eta_zeolite.value

        for w in self.W:
            total_flow = self.flow_mol_total[t, w].value
            if total_flow <= 0.0:
                continue

            pressure_w = self.pressure[t, w].value
            y = {
                comp: self.flow_mol_comp[t, w, comp].value / total_flow
                for comp in self.component_list
            }
            p = {comp: y[comp] * pressure_w / 1e5 + 1e-6 for comp in self.component_list}

            self.rate_rwgs[t, w].set_value(
                rate_mult * self.k_rwgs.value * (
                    p['CO2'] * p['H2'] - p['CO'] * p['H2O'] / self.Keq_rwgs.value
                )
            )
            if self.config.kinetics_model == 'marvast_2005':
                p_co_safe = p['CO'] + 1e-8
                p_h2_safe = p['H2'] + 1e-8
                k_c1 = self.k0_c1.value * pyo.exp(-self.E_c1.value / (8.314 * self.temperature[t, w].value))
                k_c2h4 = self.k0_c2h4.value * pyo.exp(-self.E_c2h4.value / (8.314 * self.temperature[t, w].value))
                k_c2h6 = self.k0_c2h6.value * pyo.exp(-self.E_c2h6.value / (8.314 * self.temperature[t, w].value))
                k_c3h8 = self.k0_c3h8.value * pyo.exp(-self.E_c3h8.value / (8.314 * self.temperature[t, w].value))
                k_nc4 = self.k0_nc4h10.value * pyo.exp(-self.E_nc4h10.value / (8.314 * self.temperature[t, w].value))
                k_ic4 = self.k0_ic4h10.value * pyo.exp(-self.E_ic4h10.value / (8.314 * self.temperature[t, w].value))
                k_c5 = self.k0_c5_c12.value * pyo.exp(-self.E_c5_c12.value / (8.314 * self.temperature[t, w].value))

                self.rate_c1[t, w].set_value(
                    rate_mult * k_c1 * p_co_safe**self.m_c1.value * p_h2_safe**self.n_c1.value
                )
                self.rate_c2_c4[t, w].set_value(
                    rate_mult * (
                        k_c2h4 * p_co_safe**self.m_c2h4.value * p_h2_safe**self.n_c2h4.value
                        + k_c2h6 * p_co_safe**self.m_c2h6.value * p_h2_safe**self.n_c2h6.value
                        + k_c3h8 * p_co_safe**self.m_c3h8.value * p_h2_safe**self.n_c3h8.value
                        + k_nc4 * p_co_safe**self.m_nc4h10.value * p_h2_safe**self.n_nc4h10.value
                        + k_ic4 * p_co_safe**self.m_ic4h10.value * p_h2_safe**self.n_ic4h10.value
                    )
                )
                self.rate_c5_c12[t, w].set_value(
                    rate_mult * k_c5 * p_co_safe**self.m_c5_c12.value * p_h2_safe**self.n_c5_c12.value
                )
                self.rate_c13_plus[t, w].set_value(0.0)
            elif self.config.kinetics_model == 'rwgs_2017':
                self.rate_c1[t, w].set_value(0.0)
                self.rate_c2_c4[t, w].set_value(0.0)
                self.rate_c5_c12[t, w].set_value(0.0)
                self.rate_c13_plus[t, w].set_value(0.0)
            else:
                self.rate_c1[t, w].set_value(rate_mult * self.k_c1.value * p['CO'] * p['H2'])
                self.rate_c2_c4[t, w].set_value(rate_mult * self.k_c2_c4.value * p['CO'] * p['H2'])
                self.rate_c5_c12[t, w].set_value(rate_mult * self.k_c5_c12.value * p['CO'] * p['H2'])
                self.rate_c13_plus[t, w].set_value(rate_mult * self.k_c13_plus.value * p['CO'] * p['H2'])

            if self.config.include_zeolite_reactions:
                self.rate_cracking[t, w].set_value(rate_mult * self.k_cracking.value * y['C13_plus'])
                self.rate_light_cracking[t, w].set_value(rate_mult * self.k_light_cracking.value * y['C5_C12'])
                self.rate_isomerization[t, w].set_value(rate_mult * self.k_isomerization.value * y['C5_C12'])
                self.rate_oligomerization[t, w].set_value(rate_mult * self.k_oligomerization.value * y['C2_C4'])
                self.rate_aromatization[t, w].set_value(rate_mult * self.k_aromatization.value * y['C2_C4'])
                self.rate_coke_formation[t, w].set_value(rate_mult * self.k_coke_formation.value * y['aromatics'])
            else:
                self.rate_cracking[t, w].set_value(0.0)
                self.rate_light_cracking[t, w].set_value(0.0)
                self.rate_isomerization[t, w].set_value(0.0)
                self.rate_oligomerization[t, w].set_value(0.0)
                self.rate_aromatization[t, w].set_value(0.0)
                self.rate_coke_formation[t, w].set_value(0.0)

            for comp in self.component_list:
                ft_term = sum(
                    FT_REACTIONS[r].get(comp, 0.0) * self.reaction_rate[r][t, w].value
                    for r in FT_REACTIONS
                )
                zeo_term = 0.0
                if self.config.include_zeolite_reactions:
                    zeo_term = sum(
                        ZEOLITE_REACTIONS[r].get(comp, 0.0) * self.reaction_rate[r][t, w].value
                        for r in ZEOLITE_REACTIONS
                    )
                self.dF_dW[t, w, comp].set_value(self.W_total.value * (eta_ft * ft_term + eta_zeo * zeo_term))

        print("[OK] Reactor initialized successfully")
        print(
            f"  Inlet: CO2={inlet_flow.get('CO2', 0):.3f}, H2={inlet_flow.get('H2', 0):.3f} kmol/s"
        )
        print(f"  T={temperature:.1f} K, P={pressure/1e5:.1f} bar")


def build_ft_rwgs_reactor(
    flowsheet,
    inlet_flow: Dict[str, float],
    temperature: float = 523.15,
    pressure: float = 2e6,
    W_total: float = 1.0,
    k_rwgs: float = 0.1,
    Keq_rwgs: float = 0.8,
    k_c1: float = 0.05,
    k_c2_c4: float = 0.02,
    k_c5_c12: float = 0.01,
    k_c13_plus: float = 0.005,
) -> FTRWGSReactorData:
    """
    Build and initialize an RWGS + FT reactor.
    
    Parameters
    ----------
    flowsheet : FlowsheetBlock
        Parent flowsheet
    inlet_flow : dict
        Inlet flows {component: kmol/s}
    temperature : float
        Inlet temperature [K]
    pressure : float
        Inlet pressure [Pa]
    W_total : float
        Total catalyst weight [kg]
    k_rwgs : float
        RWGS rate constant
    Keq_rwgs : float
        RWGS equilibrium constant
    k_c1 : float
        C1 formation rate constant
    k_c2_c4 : float
        C2-C4 formation rate constant
    k_c5_c12 : float
        C5-C12 formation rate constant
    k_c13_plus : float
        C13+ formation rate constant
    
    Returns
    -------
    FTRWGSReactorData
        Initialized reactor model
    """
    reactor = FTRWGSReactor()
    
    reactor.initialize(
        inlet_flow=inlet_flow,
        temperature=temperature,
        pressure=pressure,
        W_total=W_total,
        k_rwgs=k_rwgs,
        Keq_rwgs=Keq_rwgs,
        k_c1=k_c1,
        k_c2_c4=k_c2_c4,
        k_c5_c12=k_c5_c12,
        k_c13_plus=k_c13_plus,
    )
    
    return reactor


def discretize_reactor(reactor, nfe: int = 10, method: str = 'BACKWARD'):
    """
    Apply discretization to the reactor spatial domain.

    Parameters
    ----------
    reactor : FTRWGSReactorData
        The reactor to discretize
    nfe : int
        Number of finite elements
    method : str
        Discretization method ('BACKWARD', 'FORWARD', 'CENTRAL')
    """
    discretizer = TransformationFactory('dae.finite_difference')
    discretizer.apply_to(reactor, wrt=reactor.W, nfe=nfe, scheme=method)
    print(f"[OK] Reactor discretized with {nfe} finite elements ({method})")


def verify_atom_conservation(model, inlet_flow: Dict[str, float], tolerance: float = 1e-2):
    """
    Verify that atoms (C, H, O) are conserved after solving.

    Parameters
    ----------
    model : ConcreteModel
        The solved model
    inlet_flow : dict
        Inlet flows {component: kmol/s}
    tolerance : float
        Acceptable relative error for conservation (default 1%)
    """
    t = model.fs.time.first()
    reactor = model.fs.reactor

    # Inlet atoms
    inlet_atoms = count_atoms_in_stream(inlet_flow)

    # Outlet atoms
    w_outlet = max(reactor.W)
    outlet_flow = {
        c: value(reactor.flow_mol_comp[t, w_outlet, c])
        for c in reactor.component_list
    }
    outlet_atoms = count_atoms_in_stream(outlet_flow)

    print("\n" + "="*70)
    print("ATOM CONSERVATION VERIFICATION")
    print("="*70)

    print("\nINLET STREAM:")
    print(f"  Carbon:   {inlet_atoms['C']:.6f} kmol")
    print(f"  Hydrogen: {inlet_atoms['H']:.6f} kmol")
    print(f"  Oxygen:   {inlet_atoms['O']:.6f} kmol")

    print("\nOUTLET STREAM:")
    print(f"  Carbon:   {outlet_atoms['C']:.6f} kmol")
    print(f"  Hydrogen: {outlet_atoms['H']:.6f} kmol")
    print(f"  Oxygen:   {outlet_atoms['O']:.6f} kmol")

    print("\nDIFFERENCES:")
    errors = {}
    for atom in ['C', 'H', 'O']:
        diff = outlet_atoms[atom] - inlet_atoms[atom]
        rel_error = abs(diff) / (inlet_atoms[atom] + 1e-10) * 100
        errors[atom] = (diff, rel_error)
        print(f"  {atom}: {diff:+.6f} kmol ({rel_error:+.2f}%)")

    for atom, (diff, rel_error) in errors.items():
        if rel_error > tolerance:
            raise AssertionError(
                f"Atom {atom} conservation violated: {rel_error:.2f}% error > {tolerance}% tolerance"
            )

    print("\n[OK] All atoms conserved within tolerance!")
    print("="*70 + "\n")

    return True




def run_single_simulation(sim_config: Dict[str, float], inlet_flow: Dict[str, float]) -> bool:
    """
    Run a single simulation with the provided configuration.
    """
    from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition
    from idaes.core import FlowsheetBlock

    print("\n" + "="*70)
    print("RWGS + FT + ZEOLITE REACTOR - SINGLE SIMULATION")
    print("="*70)

    print("\n1. Building model...")
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0])

    print("\n2. Creating reactor...")
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

    print("\n3. Discretizing spatial domain...")
    discretize_reactor(m.fs.reactor, nfe=int(sim_config['nfe']))

    t = 0
    W_inlet = m.fs.reactor.W.first()
    for c in m.fs.reactor.component_list:
        m.fs.reactor.flow_mol_comp[t, W_inlet, c].set_value(inlet_flow.get(c, 0.0))

    m.fs.reactor.temperature[t, W_inlet].set_value(sim_config['temperature'])
    m.fs.reactor.pressure[t, W_inlet].set_value(sim_config['pressure_bar'] * 101325.0)

    print("\n4. Solving with IPOPT...")
    solver = SolverFactory('ipopt')
    solver.options['max_iter'] = int(sim_config['max_iter'])
    solver.options['tol'] = sim_config['tol']

    try:
        from pyomo.util.infeasible import log_infeasible_constraints, log_infeasible_bounds
        import logging

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
                    print(f"[WARN] Staged solve at rate_multiplier={ramp_value} failed: {stage_results.solver.termination_condition}")
                    break

            if m.fs.reactor.config.pressure_drop:
                for pd_value in (0.0, 0.1, 0.3, 0.6, 1.0):
                    m.fs.reactor.pressure_drop_multiplier.set_value(pd_value)
                    solver.options['max_iter'] = 400
                    stage_results = solver.solve(m, tee=False)
                    if stage_results.solver.termination_condition != TerminationCondition.optimal:
                        print(f"[WARN] Staged solve at pressure_drop_multiplier={pd_value} failed: {stage_results.solver.termination_condition}")
                        break

            m.fs.reactor.eta_ft.set_value(original_eta_ft)
            m.fs.reactor.eta_zeolite.set_value(original_eta_zeo)
            m.fs.reactor.ua_per_kg.set_value(original_ua)
            m.fs.reactor.rate_multiplier.set_value(original_rate_multiplier)
            m.fs.reactor.pressure_drop_multiplier.set_value(original_pressure_drop_multiplier)

            solver.options['max_iter'] = int(sim_config['max_iter'])

        results = solver.solve(m, tee=True)

        if results.solver.termination_condition == TerminationCondition.optimal:
            print("[OK] Solution converged!")

            print("\n" + "="*70)
            print("REACTOR RESULTS")
            print("="*70)

            print("\nINLET (W=0):")
            for c in m.fs.reactor.component_list:
                flow = value(m.fs.reactor.flow_mol_comp[t, W_inlet, c])
                print(f"  {c:8s}: {flow:.6f} kmol/s")

            print("\nOUTLET (W=1):")
            W_outlet = m.fs.reactor.W.last()
            for c in m.fs.reactor.component_list:
                flow = value(m.fs.reactor.flow_mol_comp[t, W_outlet, c])
                print(f"  {c:8s}: {flow:.6f} kmol/s")

            inlet_CO = value(m.fs.reactor.flow_mol_comp[t, W_inlet, 'CO'])
            outlet_CO = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'CO'])
            CO_conversion = 100 * (inlet_CO - outlet_CO) / inlet_CO if inlet_CO > 1e-6 else 0

            inlet_CO2 = value(m.fs.reactor.flow_mol_comp[t, W_inlet, 'CO2'])
            outlet_CO2 = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'CO2'])
            CO2_conversion = 100 * (inlet_CO2 - outlet_CO2) / inlet_CO2 if inlet_CO2 > 1e-6 else 0

            inlet_H2 = value(m.fs.reactor.flow_mol_comp[t, W_inlet, 'H2'])
            outlet_H2 = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'H2'])
            H2_conversion = 100 * (inlet_H2 - outlet_H2) / inlet_H2 if inlet_H2 > 1e-6 else 0

            inlet_H2O = value(m.fs.reactor.flow_mol_comp[t, W_inlet, 'H2O'])
            outlet_H2O = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'H2O'])
            H2O_formation = outlet_H2O - inlet_H2O

            print(f"\nCO2 Conversion: {CO2_conversion:.2f}%")
            print(f"\nCO Conversion: {CO_conversion:.2f}%")
            print(f"H2 Conversion: {H2_conversion:.2f}%")
            print(f"H2O Formation: {H2O_formation:.6f} kmol/s")

            print("\nProduct Formation:")
            c1_out = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'C1'])
            c2c4_out = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'C2_C4'])
            c5c12_out = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'C5_C12'])
            c13p_out = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'C13_plus'])
            print(f"  C1 yield:        {c1_out:.6f} kmol/s")
            print(f"  C2-C4 yield:     {c2c4_out:.6f} kmol/s")
            print(f"  C5-C12 yield:    {c5c12_out:.6f} kmol/s")
            print(f"  C13+ yield:      {c13p_out:.6f} kmol/s")
            print("="*70)

            print("\n5. Verifying atom conservation...")
            verify_atom_conservation(m, inlet_flow, tolerance=1.0)

            if m.fs.reactor.config.energy_balance:
                print("\nTemperature profile (K):")
                for w in m.fs.reactor.W:
                    T_w = value(m.fs.reactor.temperature[t, w])
                    print(f"  W={w:.3f}: T={T_w:.2f}")

            return True

        print(f"[FAIL] Solver failed: {results.solver.termination_condition}")
        logging.getLogger('pyomo.util.infeasible').setLevel(logging.INFO)
        log_infeasible_constraints(m, tol=1e-6)
        log_infeasible_bounds(m, tol=1e-6)
        return False
    except Exception as e:
        print(f"[FAIL] Error during solve: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    # ==================== EDIT PARAMETERS HERE ====================
    from sim_config import SIM_CONFIG

    # Inlet composition (kmol/s): CO2 + H2 feed
    INLET_FLOW = {
        'CO2': 0.3,
        'H2': 0.7,
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
    # ==============================================================

    run_single_simulation(SIM_CONFIG, INLET_FLOW)

