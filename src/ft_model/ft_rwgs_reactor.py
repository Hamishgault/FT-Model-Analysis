"""
RWGS + Fischer-Tropsch Packed-Bed Reactor

This reactor combines the water-gas shift reaction with Fischer-Tropsch
product synthesis in a single catalytic bed.

Reactions:
1. RWGS:    CO2 + H2   <-> CO   + H2O
2. CH4:     CO  + 2H2  <-> CH4  + H2O
3. C2H4:    CO  + 2H2  <-> C2H4 + H2O
4. C5+:     CO  + 2H2  <-> C5+  + H2O

This version:
- Uses consistent indexing with expanded component list
- All reactions use power-law rate expressions
- Can be discretized and solved with IPOPT
- Demonstrates multi-reaction system with proper stoichiometry
"""
# type: ignore  # Pyomo/IDAES type hints not fully recognized by Pylance

from typing import Dict, Optional

import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.environ import TransformationFactory, value

from idaes.core import (
    UnitModelBlockData,
    declare_process_block_class,
)


# ==================== STOICHIOMETRIC COEFFICIENTS ====================

# RWGS: CO2 + H2 <-> CO + H2O
RWGS_STOICHIOMETRY = {
    'CO2': -1.0,
    'H2': -1.0,
    'CO': 1.0,
    'H2O': 1.0,
    'CH4': 0.0,
    'C2H4': 0.0,
    'C5plus': 0.0,
}

# Methane formation: CO + 2H2 -> CH4 + H2O
CH4_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -2.0,
    'CO': -1.0,
    'H2O': 1.0,
    'CH4': 1.0,
    'C2H4': 0.0,
    'C5plus': 0.0,
}

# Ethylene formation: CO + 2H2 -> C2H4 + H2O
C2H4_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -2.0,
    'CO': -1.0,
    'H2O': 1.0,
    'CH4': 0.0,
    'C2H4': 1.0,
    'C5plus': 0.0,
}

# C5+ lump formation: CO + 2H2 -> C5+ + H2O
C5PLUS_STOICHIOMETRY = {
    'CO2': 0.0,
    'H2': -2.0,
    'CO': -1.0,
    'H2O': 1.0,
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


@declare_process_block_class('FTRWGSReactor')
class FTRWGSReactorData(UnitModelBlockData):
    """
    RWGS + Fischer-Tropsch packed-bed reactor with proper indexing and DAE solving.
    
    Reactions:
    1. RWGS:    CO2 + H2   <-> CO   + H2O
    2. CH4:     CO  + 2H2  <-> CH4  + H2O
    3. C2H4:    CO  + 2H2  <-> C2H4 + H2O
    4. C5+:     CO  + 2H2  <-> C5+  + H2O
    
    This reactor uses 1D spatial discretization along catalyst weight (W).
    All indexing is consistent and the model can be discretized and solved.
    
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
        
        # CH4 formation rate constant (power-law)
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
        
        # CH4 formation rate: r = k * p_CO * p_H2^2 (power-law)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='CH4 formation rate',
        )
        def rate_ch4_eq(b, t, w):
            p_CO = b.partial_pressure[t, w, 'CO']
            p_H2 = b.partial_pressure[t, w, 'H2']
            return b.rate_ch4[t, w] == b.k_ch4 * p_CO * (p_H2 ** 2)
        
        # C2H4 formation rate: r = k * p_CO * p_H2^2
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='C2H4 formation rate',
        )
        def rate_c2h4_eq(b, t, w):
            p_CO = b.partial_pressure[t, w, 'CO']
            p_H2 = b.partial_pressure[t, w, 'H2']
            return b.rate_c2h4[t, w] == b.k_c2h4 * p_CO * (p_H2 ** 2)
        
        # C5+ formation rate: r = k * p_CO * p_H2^2
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='C5+ formation rate',
        )
        def rate_c5plus_eq(b, t, w):
            p_CO = b.partial_pressure[t, w, 'CO']
            p_H2 = b.partial_pressure[t, w, 'H2']
            return b.rate_c5plus[t, w] == b.k_c5plus * p_CO * (p_H2 ** 2)
        
        # Material balance: dF_i/dW = sum over reactions of (nu_ij * r_j)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            self.component_list,
            doc='Material balance ODEs',
        )
        def material_balance(b, t, w, c):
            if w == 0:
                return pyo.Constraint.Skip  # Boundary condition, not ODE
            
            # Sum stoichiometric contributions from all reactions
            rhs = (
                RWGS_STOICHIOMETRY[c] * b.rate_rwgs[t, w] +
                CH4_STOICHIOMETRY[c] * b.rate_ch4[t, w] +
                C2H4_STOICHIOMETRY[c] * b.rate_c2h4[t, w] +
                C5PLUS_STOICHIOMETRY[c] * b.rate_c5plus[t, w]
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
        
        print("[OK] Reactor initialized successfully")
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


if __name__ == '__main__':
    from pyomo.environ import ConcreteModel, SolverFactory
    from idaes.core import FlowsheetBlock
    
    print("\n" + "="*70)
    print("RWGS + FT REACTOR - SELF TEST")
    print("="*70)
    
    # Create model
    print("\n1. Building model...")
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0])
    
    # Define inlet conditions: H2 and CO synthesis gas (product from RWGS)
    inlet_flow = {
        'CO2': 1e-8,     # Trace amounts to satisfy bounds
        'H2': 0.4,
        'CO': 0.6,
        'H2O': 1e-8,
        'CH4': 1e-8,
        'C2H4': 1e-8,
        'C5plus': 1e-8,
    }
    
    # Build reactor
    print("\n2. Creating reactor...")
    m.fs.reactor = FTRWGSReactor()
    m.fs.reactor.initialize(
        inlet_flow=inlet_flow,
        temperature=523.15,  # 250°C
        pressure=20e5,  # 20 bar
        W_total=5.0,  # 5 kg catalyst
        k_rwgs=0.05,
        Keq_rwgs=0.8,
        k_ch4=0.02,      # Reduce CH4 rate for stability
        k_c2h4=0.01,     # Reduce C2H4 rate
        k_c5plus=0.005,  # Reduce C5+ rate
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
            print("RESULTS")
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
        else:
            print(f"[FAIL] Solver failed: {results.solver.termination_condition}")
    except Exception as e:
        print(f"[FAIL] Error during solve: {e}")
        import traceback
        traceback.print_exc()
