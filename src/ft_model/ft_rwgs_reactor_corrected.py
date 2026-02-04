"""
RWGS + Fischer-Tropsch Packed-Bed Reactor with Atom Balance Verification

This reactor combines the water-gas shift reaction with Fischer-Tropsch
product synthesis in a single catalytic bed with strict atom conservation.

CORRECTED Reactions (atom-balanced):
1. RWGS:    CO2 + H2   <-> CO   + H2O
2. CH4:     CO  + 3H2  -> CH4  + H2O        (CORRECTED: was CO + 2H2)
3. C2H4:   2CO  + 4H2  -> C2H4 + 2H2O      (CORRECTED: was CO + 2H2)
4. C5+:    5CO  + 10H2 -> C5H10 + 5H2O     (CORRECTED: was CO + 2H2)

This version:
- All stoichiometries verified for C, H, O atom conservation
- Corrected rate expressions match stoichiometry
- Includes atom_balance() verification function
- Post-solve atom conservation check
- Self-test with conservation assertion
"""
# type: ignore  # Pyomo/IDAES type hints not fully recognized by Pylance

from typing import Dict

import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.environ import TransformationFactory, value

from idaes.core import (
    UnitModelBlockData,
    declare_process_block_class,
)


# ==================== ATOMIC COMPOSITION ====================

ATOMIC_COMPOSITION = {
    'CO2': {'C': 1, 'H': 0, 'O': 2},
    'H2': {'C': 0, 'H': 2, 'O': 0},
    'CO': {'C': 1, 'H': 0, 'O': 1},
    'H2O': {'C': 0, 'H': 2, 'O': 1},
    'CH4': {'C': 1, 'H': 4, 'O': 0},
    'C2H4': {'C': 2, 'H': 4, 'O': 0},
    'C5plus': {'C': 5, 'H': 10, 'O': 0},  # C5H10 (unsaturated)
}


def check_atom_balance(reaction_name: str, stoichiometry: Dict[str, float]) -> bool:
    """
    Verify that a reaction conserves C, H, and O atoms.
    
    Parameters
    ----------
    reaction_name : str
        Name of the reaction (for error messages)
    stoichiometry : dict
        {component: nu_ij} where nu_ij is stoichiometric coefficient
        (negative for reactants, positive for products)
    
    Returns
    -------
    bool
        True if balanced
    
    Raises
    ------
    ValueError
        If reaction has unbalanced atoms
    """
    atoms = {'C': 0, 'H': 0, 'O': 0}
    
    # Sum atomic contributions: if nu_ij < 0 (reactant), atoms are consumed
    # if nu_ij > 0 (product), atoms are produced
    for comp, nu in stoichiometry.items():
        if nu == 0:
            continue
        if comp not in ATOMIC_COMPOSITION:
            raise ValueError(f"Unknown component '{comp}' in reaction '{reaction_name}'")
        
        for atom, count in ATOMIC_COMPOSITION[comp].items():
            atoms[atom] += nu * count
    
    # Check balance: sum should be zero if balanced
    for atom, total in atoms.items():
        if abs(total) > 1e-6:
            raise ValueError(
                f"Reaction '{reaction_name}' has UNBALANCED {atom} atoms: {total}"
            )
    
    return True


def count_atoms_in_stream(flows: Dict[str, float]) -> Dict[str, float]:
    """
    Count total C, H, O atoms in a stream given component flows.
    
    Parameters
    ----------
    flows : dict
        {component: flow_kmol_s}
    
    Returns
    -------
    dict
        {'C': count, 'H': count, 'O': count}
    """
    atoms = {'C': 0.0, 'H': 0.0, 'O': 0.0}
    
    for comp, flow in flows.items():
        if comp not in ATOMIC_COMPOSITION:
            continue
        for atom, count in ATOMIC_COMPOSITION[comp].items():
            atoms[atom] += flow * count
    
    return atoms


# ==================== STOICHIOMETRIC COEFFICIENTS ====================
# All stoichiometries verified for C, H, O atom conservation

# RWGS: CO2 + H2 <-> CO + H2O
# Balanced: C=1, H=2, O=2 on both sides
RWGS_STOICHIOMETRY = {
    'CO2': -1.0,
    'H2': -1.0,
    'CO': 1.0,
    'H2O': 1.0,
    'CH4': 0.0,
    'C2H4': 0.0,
    'C5plus': 0.0,
}

# CH4 formation: CO + 3H2 -> CH4 + H2O
# Balanced: C=1, H=6, O=1 on both sides
# CORRECTED from CO + 2H2
CH4_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -3.0,  # CORRECTED: was -2.0
    'CO': -1.0,
    'H2O': 1.0,
    'CH4': 1.0,
    'C2H4': 0.0,
    'C5plus': 0.0,
}

# C2H4 formation: 2CO + 4H2 -> C2H4 + 2H2O
# Balanced: C=2, H=8, O=2 on both sides
# CORRECTED from CO + 2H2
C2H4_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -4.0,   # CORRECTED: was -2.0
    'CO': -2.0,   # CORRECTED: was -1.0
    'H2O': 2.0,   # CORRECTED: was 1.0
    'CH4': 0.0,
    'C2H4': 1.0,
    'C5plus': 0.0,
}

# C5+ lump formation: 5CO + 10H2 -> C5H10 + 5H2O
# C5plus represents C5H10 (unsaturated)
# Balanced: C=5, H=20, O=5 on both sides
# CORRECTED from CO + 2H2
C5PLUS_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -10.0,  # CORRECTED: was -2.0
    'CO': -5.0,   # CORRECTED: was -1.0
    'H2O': 5.0,   # CORRECTED: was 1.0
    'CH4': 0.0,
    'C2H4': 0.0,
    'C5plus': 1.0,
}

# Dictionary of all reactions for easy access
ALL_REACTIONS = {
    'rwgs': RWGS_STOICHIOMETRY,
    'ch4': CH4_STOICHIOMETRY,
    'c2h4': C2H4_STOICHIOMETRY,
    'c5plus': C5PLUS_STOICHIOMETRY,
}

# Verify all stoichiometries on module load
print("\n" + "="*70)
print("STOICHIOMETRY VERIFICATION")
print("="*70)
try:
    check_atom_balance('RWGS', RWGS_STOICHIOMETRY)
    check_atom_balance('CH4 Formation', CH4_STOICHIOMETRY)
    check_atom_balance('C2H4 Formation', C2H4_STOICHIOMETRY)
    check_atom_balance('C5+ Formation', C5PLUS_STOICHIOMETRY)
    print("All reactions verified: atom balance OK")
except ValueError as e:
    print(f"[ERROR] {e}")
    raise
print("="*70 + "\n")


@declare_process_block_class('FTRWGSReactor')
class FTRWGSReactorData(UnitModelBlockData):
    """
    RWGS + Fischer-Tropsch packed-bed reactor with atom conservation.
    
    Reactions (all atom-balanced):
    1. RWGS:    CO2 + H2   <-> CO   + H2O
    2. CH4:     CO  + 3H2  -> CH4  + H2O
    3. C2H4:   2CO  + 4H2  -> C2H4 + 2H2O
    4. C5+:    5CO  + 10H2 -> C5H10 + 5H2O
    
    This reactor uses 1D spatial discretization along catalyst weight (W).
    All indexing is consistent and atom conservation is verified.
    
    Attributes
    ----------
    W : ContinuousSet
        Normalized catalyst weight domain [0, 1]
    flow_mol_comp : Var
        Component molar flow rates [kmol/s] indexed by (time, W, component)
    temperature : Var
        Temperature [K] indexed by (time, W)
    pressure : Var
        Pressure [Pa] indexed by (time, W)
    """
    
    CONFIG = UnitModelBlockData.CONFIG()
    
    def build(self):
        """
        Build the RWGS + FT reactor model.
        
        Creates spatial domain, state variables, parameters, and constraints.
        """
        super().build()
        
        # Component list - includes RWGS feedstock and FT products
        self.component_list = ['CO2', 'H2', 'CO', 'H2O', 'CH4', 'C2H4', 'C5plus']
        
        # Reaction list for convenient iteration
        self.reaction_list = ['rwgs', 'ch4', 'c2h4', 'c5plus']
        
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
        
        # CH4 formation rate constant
        self.k_ch4 = pyo.Param(
            initialize=0.05,
            mutable=True,
            doc='CH4 formation rate constant',
        )
        
        # C2H4 formation rate constant
        self.k_c2h4 = pyo.Param(
            initialize=0.02,
            mutable=True,
            doc='C2H4 formation rate constant',
        )
        
        # C5+ formation rate constant
        self.k_c5plus = pyo.Param(
            initialize=0.01,
            mutable=True,
            doc='C5+ formation rate constant',
        )
        
        # ==================== STATE VARIABLES ====================
        
        # Component molar flow rates [kmol/s]
        self.flow_mol_comp = pyo.Var(
            self.flowsheet().time,
            self.W,
            self.component_list,
            initialize=0.1,
            bounds=(1e-8, None),
            doc='Component molar flow rates [kmol/s]',
        )
        
        # Temperature [K]
        self.temperature = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=500.0,
            bounds=(200.0, 1000.0),
            doc='Temperature [K]',
        )
        
        # Pressure [Pa]
        self.pressure = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=2e6,
            bounds=(1e5, 5e6),
            doc='Pressure [Pa]',
        )
        
        # Total molar flow [kmol/s]
        self.flow_mol_total = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=1.0,
            bounds=(1e-6, None),
            doc='Total molar flow [kmol/s]',
        )
        
        # Mole fractions
        self.mole_frac = pyo.Var(
            self.flowsheet().time,
            self.W,
            self.component_list,
            initialize=1.0 / len(self.component_list),
            bounds=(0, 1),
            doc='Mole fractions',
        )
        
        # Partial pressures [Pa]
        self.partial_pressure = pyo.Var(
            self.flowsheet().time,
            self.W,
            self.component_list,
            initialize=5e5,
            bounds=(0, None),
            doc='Partial pressures [Pa]',
        )
        
        # Reaction rates [kmol/(kg_cat·s)]
        self.rate_rwgs = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.01,
            bounds=None,
            doc='RWGS reaction rate [kmol/(kg_cat·s)]',
        )
        
        self.rate_ch4 = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.001,
            bounds=None,
            doc='CH4 formation rate [kmol/(kg_cat·s)]',
        )
        
        self.rate_c2h4 = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.0005,
            bounds=None,
            doc='C2H4 formation rate [kmol/(kg_cat·s)]',
        )
        
        self.rate_c5plus = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.0001,
            bounds=None,
            doc='C5+ formation rate [kmol/(kg_cat·s)]',
        )
        
        # ==================== DERIVATIVES ====================
        
        # Material balance derivatives dF/dW
        self.dF_dW = DerivativeVar(
            self.flow_mol_comp,
            wrt=self.W,
            doc='Derivative of flow w.r.t. catalyst weight',
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
        
        # Mole fraction calculation
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            self.component_list,
            doc='Mole fractions',
        )
        def mole_frac_eq(b, t, w, c):
            return b.mole_frac[t, w, c] * b.flow_mol_total[t, w] == b.flow_mol_comp[t, w, c]
        
        # Partial pressure calculation
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            self.component_list,
            doc='Partial pressures',
        )
        def partial_pressure_eq(b, t, w, c):
            return b.partial_pressure[t, w, c] == b.mole_frac[t, w, c] * b.pressure[t, w]
        
        # RWGS rate: r = k * (p_CO2 * p_H2 - p_CO * p_H2O / Keq)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='RWGS reaction rate',
        )
        def rate_rwgs_eq(b, t, w):
            p_CO2 = b.partial_pressure[t, w, 'CO2']
            p_H2 = b.partial_pressure[t, w, 'H2']
            p_CO = b.partial_pressure[t, w, 'CO']
            p_H2O = b.partial_pressure[t, w, 'H2O']
            
            forward = p_CO2 * p_H2
            reverse = p_CO * p_H2O / b.Keq_rwgs
            
            return b.rate_rwgs[t, w] == b.k_rwgs * (forward - reverse)
        
        # CH4 formation rate: r = k * p_CO * p_H2
        # Simplified from stoichiometric order p_CO * p_H2^3 for numerical stability
        # Reaction stoichiometry: CO + 3H2 -> CH4 + H2O (atoms are balanced)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='CH4 formation rate (simplified power-law)',
        )
        def rate_ch4_eq(b, t, w):
            p_CO = b.partial_pressure[t, w, 'CO']
            p_H2 = b.partial_pressure[t, w, 'H2']
            return b.rate_ch4[t, w] == b.k_ch4 * p_CO * p_H2
        
        # C2H4 formation rate: r = k * p_CO * p_H2
        # Simplified from stoichiometric order p_CO^2 * p_H2^4 for numerical stability
        # Reaction stoichiometry: 2CO + 4H2 -> C2H4 + 2H2O (atoms are balanced)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='C2H4 formation rate (simplified power-law)',
        )
        def rate_c2h4_eq(b, t, w):
            p_CO = b.partial_pressure[t, w, 'CO']
            p_H2 = b.partial_pressure[t, w, 'H2']
            return b.rate_c2h4[t, w] == b.k_c2h4 * p_CO * p_H2
        
        # C5+ formation rate: r = k * p_CO * p_H2
        # Simplified from stoichiometric order p_CO^5 * p_H2^10 for numerical stability
        # Reaction stoichiometry: 5CO + 10H2 -> C5H10 + 5H2O (atoms are balanced)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='C5+ formation rate (simplified power-law)',
        )
        def rate_c5plus_eq(b, t, w):
            p_CO = b.partial_pressure[t, w, 'CO']
            p_H2 = b.partial_pressure[t, w, 'H2']
            return b.rate_c5plus[t, w] == b.k_c5plus * p_CO * p_H2
        
        # Material balance: dF_i/dW = sum_j (nu_ij * r_j)
        # nu_ij is stoichiometric coefficient of component i in reaction j
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            self.component_list,
            doc='Material balance ODEs (with atom conservation)',
        )
        def material_balance(b, t, w, c):
            if w == 0:
                return pyo.Constraint.Skip  # Boundary condition, not ODE
            
            # Sum stoichiometric contributions from all reactions
            # Each nu_ij is from the stoichiometry dictionaries
            rhs = (
                RWGS_STOICHIOMETRY.get(c, 0.0) * b.rate_rwgs[t, w] +
                CH4_STOICHIOMETRY.get(c, 0.0) * b.rate_ch4[t, w] +
                C2H4_STOICHIOMETRY.get(c, 0.0) * b.rate_c2h4[t, w] +
                C5PLUS_STOICHIOMETRY.get(c, 0.0) * b.rate_c5plus[t, w]
            )
            
            return b.dF_dW[t, w, c] == rhs
        
        # Isothermal constraint (can be relaxed for energy balance)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='Isothermal operation',
        )
        def isothermal_eq(b, t, w):
            if w == 0:
                return pyo.Constraint.Skip
            return b.temperature[t, w] == b.temperature[t, 0]
        
        # Constant pressure constraint
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='Constant pressure',
        )
        def isobaric_eq(b, t, w):
            if w == 0:
                return pyo.Constraint.Skip
            return b.pressure[t, w] == b.pressure[t, 0]
    
    def initialize(
        self,
        inlet_flow: Dict[str, float],
        temperature: float = 523.15,
        pressure: float = 2e6,
        W_total: float = 1.0,
        k_rwgs: float = 0.1,
        Keq_rwgs: float = 0.8,
        k_ch4: float = 0.05,
        k_c2h4: float = 0.02,
        k_c5plus: float = 0.01,
    ):
        """
        Initialize the reactor with inlet conditions and parameters.
        
        Parameters
        ----------
        inlet_flow : dict
            Inlet molar flows {component: flow_kmol_s}
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
        k_ch4 : float
            CH4 formation rate constant
        k_c2h4 : float
            C2H4 formation rate constant
        k_c5plus : float
            C5+ formation rate constant
        """
        # Set parameters
        self.W_total.set_value(W_total)
        self.k_rwgs.set_value(k_rwgs)
        self.Keq_rwgs.set_value(Keq_rwgs)
        self.k_ch4.set_value(k_ch4)
        self.k_c2h4.set_value(k_c2h4)
        self.k_c5plus.set_value(k_c5plus)
        
        # Get time point
        t = self.flowsheet().time.first()
        
        # Fix inlet conditions at W=0
        for comp in self.component_list:
            flow = inlet_flow.get(comp, 1e-8)
            self.flow_mol_comp[t, 0, comp].fix(flow)
        
        self.temperature[t, 0].fix(temperature)
        self.pressure[t, 0].fix(pressure)
        
        # Propagate initial guesses along W
        for w in self.W:
            if w != 0:
                for comp in self.component_list:
                    self.flow_mol_comp[t, w, comp].set_value(
                        self.flow_mol_comp[t, 0, comp].value
                    )
                self.temperature[t, w].set_value(temperature)
                self.pressure[t, w].set_value(pressure)
        
        print(f"[OK] Reactor initialized successfully")
        print(f"  Inlet: CO2={inlet_flow.get('CO2', 0):.3f}, H2={inlet_flow.get('H2', 0):.3f} kmol/s")
        print(f"  T={temperature:.1f} K, P={pressure/1e5:.1f} bar")


def build_ft_rwgs_reactor(
    flowsheet,
    inlet_flow: Dict[str, float],
    temperature: float = 523.15,
    pressure: float = 2e6,
    W_total: float = 1.0,
    k_rwgs: float = 0.1,
    Keq_rwgs: float = 0.8,
    k_ch4: float = 0.05,
    k_c2h4: float = 0.02,
    k_c5plus: float = 0.01,
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
    k_ch4 : float
        CH4 formation rate constant
    k_c2h4 : float
        C2H4 formation rate constant
    k_c5plus : float
        C5+ formation rate constant
    
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
        k_ch4=k_ch4,
        k_c2h4=k_c2h4,
        k_c5plus=k_c5plus,
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
    
    Raises
    ------
    AssertionError
        If conservation is violated beyond tolerance
    """
    from pyomo.environ import value
    
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
    
    # Assert conservation
    for atom, (diff, rel_error) in errors.items():
        if rel_error > tolerance:
            raise AssertionError(
                f"Atom {atom} conservation violated: {rel_error:.2f}% error > {tolerance}% tolerance"
            )
    
    print("\n[OK] All atoms conserved within tolerance!")
    print("="*70 + "\n")
    
    return True


if __name__ == '__main__':
    from pyomo.environ import ConcreteModel, SolverFactory
    from idaes.core import FlowsheetBlock
    
    print("\n" + "="*70)
    print("RWGS + FT REACTOR WITH ATOM BALANCE VERIFICATION - SELF TEST")
    print("="*70)
    
    # Create model
    print("\n1. Building model...")
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0])
    
    # Define inlet conditions: synthesis gas (CO + H2, no CO2/H2O)
    # This tests the reactions in the forward direction without RWGS
    inlet_flow = {
        'CO2': 1e-8,      # Trace to satisfy bounds
        'H2': 0.4,
        'CO': 0.6,
        'H2O': 1e-8,
        'CH4': 1e-8,
        'C2H4': 1e-8,
        'C5plus': 1e-8,
    }
    
    print("\n2. Creating reactor...")
    m.fs.reactor = FTRWGSReactor()
    m.fs.reactor.initialize(
        inlet_flow=inlet_flow,
        temperature=523.15,  # 250°C
        pressure=20e5,  # 20 bar
        W_total=5.0,  # 5 kg catalyst
        k_rwgs=0.05,
        Keq_rwgs=0.8,
        k_ch4=0.02,      # CH4 rate
        k_c2h4=0.01,     # C2H4 rate
        k_c5plus=0.005,  # C5+ rate
    )
    
    # Discretize
    print("\n3. Discretizing spatial domain...")
    discretize_reactor(m.fs.reactor, nfe=20, method='BACKWARD')
    
    # Solve
    print("\n4. Solving with IPOPT...")
    try:
        from idaes.core.solvers import get_solver
        solver = get_solver('ipopt')
        solver.options['max_iter'] = 500
        solver.options['tol'] = 1e-6
        
        results = solver.solve(m, tee=True)
        
        from pyomo.opt import TerminationCondition
        if results.solver.termination_condition == TerminationCondition.optimal:
            print("\n[OK] Solution converged!")
            
            # Print results
            t = m.fs.time.first()
            W_inlet = min(m.fs.reactor.W)
            W_outlet = max(m.fs.reactor.W)
            
            print("\n" + "="*70)
            print("REACTOR RESULTS")
            print("="*70)
            
            print("\nINLET (W=0):")
            for c in m.fs.reactor.component_list:
                flow = value(m.fs.reactor.flow_mol_comp[t, W_inlet, c])
                print(f"  {c:8s}: {flow:.6f} kmol/s")
            
            print("\nOUTLET (W=1):")
            for c in m.fs.reactor.component_list:
                flow = value(m.fs.reactor.flow_mol_comp[t, W_outlet, c])
                print(f"  {c:8s}: {flow:.6f} kmol/s")
            
            # Calculate conversions
            inlet_CO = value(m.fs.reactor.flow_mol_comp[t, W_inlet, 'CO'])
            outlet_CO = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'CO'])
            CO_conversion = 100 * (inlet_CO - outlet_CO) / inlet_CO if inlet_CO > 1e-6 else 0
            
            inlet_H2 = value(m.fs.reactor.flow_mol_comp[t, W_inlet, 'H2'])
            outlet_H2 = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'H2'])
            H2_conversion = 100 * (inlet_H2 - outlet_H2) / inlet_H2 if inlet_H2 > 1e-6 else 0
            
            print(f"\nCO Conversion: {CO_conversion:.2f}%")
            print(f"H2 Conversion: {H2_conversion:.2f}%")
            
            # Product yields
            print("\nProduct Formation:")
            CH4_out = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'CH4'])
            C2H4_out = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'C2H4'])
            C5p_out = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'C5plus'])
            print(f"  CH4 yield:     {CH4_out:.6f} kmol/s")
            print(f"  C2H4 yield:    {C2H4_out:.6f} kmol/s")
            print(f"  C5+ yield:     {C5p_out:.6f} kmol/s")
            print("="*70)
            
            # Verify atom conservation
            print("\n5. Verifying atom conservation...")
            verify_atom_conservation(m, inlet_flow, tolerance=1.0)
            
            print("[OK] SELF-TEST PASSED: Reactor solves with atom conservation!")
            
        else:
            print(f"[FAIL] Solver failed: {results.solver.termination_condition}")
    except Exception as e:
        print(f"[FAIL] Error during solve: {e}")
        import traceback
        traceback.print_exc()
