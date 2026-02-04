"""
Simplified RWGS Packed-Bed Reactor - Working Implementation

This is a simplified version of the bifunctional reactor that focuses ONLY on
the RWGS reaction to ensure proper building and solving.

Reaction: CO2 + H2 <-> CO + H2O

This version:
- Uses consistent indexing
- Builds without errors
- Can be discretized and solved with IPOPT
- Serves as a template for adding FT and zeolite reactions later
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


# Stoichiometric coefficients for RWGS reaction: CO2 + H2 <-> CO + H2O
RWGS_STOICHIOMETRY = {
    'CO2': -1.0,
    'H2': -1.0,
    'CO': 1.0,
    'H2O': 1.0,
}


@declare_process_block_class('SimplifiedRWGSReactor')
class SimplifiedRWGSReactorData(UnitModelBlockData):
    """
    Simplified RWGS packed-bed reactor with proper indexing and DAE solving.
    
    Reaction: CO2 + H2 <-> CO + H2O
    
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
        Build the RWGS reactor model.
        
        Creates spatial domain, state variables, parameters, and constraints.
        """
        super().build()
        
        # Component list - RWGS only
        self.component_list = ['CO2', 'H2', 'CO', 'H2O']
        
        # Spatial domain: normalized catalyst weight [0, 1]
        self.W = ContinuousSet(bounds=(0, 1.0))
        
        # ==================== PARAMETERS ====================
        self.W_total = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Total catalyst weight [kg]',
        )
        
        self.k_rwgs = pyo.Param(
            initialize=0.1,
            mutable=True,
            doc='RWGS forward rate constant [kmol/(kg_cat·s·Pa^2)]',
        )
        
        self.Keq_rwgs = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='RWGS equilibrium constant (dimensionless)',
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
        
        # Temperature [K] - isothermal for simplicity
        self.temperature = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=500.0,
            bounds=(200.0, 1000.0),
            doc='Temperature [K]',
        )
        
        # Pressure [Pa] - constant pressure for simplicity
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
        
        # Mole fractions (dimensionless)
        self.mole_frac = pyo.Var(
            self.flowsheet().time,
            self.W,
            self.component_list,
            initialize=0.25,
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
        
        # RWGS reaction rate [kmol/(kg_cat·s)]
        self.rate_rwgs = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.01,
            bounds=None,
            doc='RWGS reaction rate [kmol/(kg_cat·s)]',
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
            doc='Total molar flow'
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
            doc='Mole fractions'
        )
        def mole_frac_eq(b, t, w, c):
            return b.mole_frac[t, w, c] * b.flow_mol_total[t, w] == b.flow_mol_comp[t, w, c]
        
        # Partial pressure calculation
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            self.component_list,
            doc='Partial pressures'
        )
        def partial_pressure_eq(b, t, w, c):
            return b.partial_pressure[t, w, c] == b.mole_frac[t, w, c] * b.pressure[t, w]
        
        # RWGS rate expression: r = k * (p_CO2 * p_H2 - p_CO * p_H2O / Keq)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='RWGS reaction rate'
        )
        def rate_rwgs_eq(b, t, w):
            p_CO2 = b.partial_pressure[t, w, 'CO2']
            p_H2 = b.partial_pressure[t, w, 'H2']
            p_CO = b.partial_pressure[t, w, 'CO']
            p_H2O = b.partial_pressure[t, w, 'H2O']
            
            forward = p_CO2 * p_H2
            reverse = p_CO * p_H2O / b.Keq_rwgs
            
            return b.rate_rwgs[t, w] == b.k_rwgs * b.W_total * (forward - reverse)
        
        # Material balance: dF_i/dW = W_total * nu_i * r_rwgs
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            self.component_list,
            doc='Material balance ODEs'
        )
        def material_balance(b, t, w, c):
            if w == 0:
                return pyo.Constraint.Skip  # Boundary condition, not ODE
            
            nu = RWGS_STOICHIOMETRY.get(c, 0.0)
            return b.dF_dW[t, w, c] == nu * b.rate_rwgs[t, w]
        
        # Isothermal constraint (can be relaxed later for energy balance)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='Isothermal operation'
        )
        def isothermal_eq(b, t, w):
            if w == 0:
                return pyo.Constraint.Skip
            return b.temperature[t, w] == b.temperature[t, 0]
        
        # Constant pressure constraint (can be relaxed for pressure drop)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            doc='Constant pressure'
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
        Keq_rwgs: float = 1.0,
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
        """
        # Set parameters
        self.W_total.set_value(W_total)
        self.k_rwgs.set_value(k_rwgs)
        self.Keq_rwgs.set_value(Keq_rwgs)
        
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
            if w > 0:
                self.temperature[t, w].set_value(temperature)
                self.pressure[t, w].set_value(pressure)
                
                for comp in self.component_list:
                    flow = inlet_flow.get(comp, 1e-8)
                    self.flow_mol_comp[t, w, comp].set_value(flow)
        
        print("✓ Reactor initialized successfully")
        print(f"  Inlet: CO2={inlet_flow.get('CO2', 0):.3f}, H2={inlet_flow.get('H2', 0):.3f} kmol/s")
        print(f"  T={temperature:.1f} K, P={pressure/1e5:.1f} bar")


def build_simplified_rwgs_reactor(
    flowsheet,
    inlet_flow: Dict[str, float],
    temperature: float = 523.15,
    pressure: float = 2e6,
    W_total: float = 1.0,
    k_rwgs: float = 0.1,
    Keq_rwgs: float = 1.0,
) -> SimplifiedRWGSReactorData:
    """
    Build and initialize a simplified RWGS reactor.
    
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
    
    Returns
    -------
    SimplifiedRWGSReactorData
        Initialized reactor model
    """
    # Create reactor - IDAES process blocks are added as attributes
    reactor = SimplifiedRWGSReactor()
    
    # Initialize
    reactor.initialize(
        inlet_flow=inlet_flow,
        temperature=temperature,
        pressure=pressure,
        W_total=W_total,
        k_rwgs=k_rwgs,
        Keq_rwgs=Keq_rwgs,
    )
    
    return reactor


def discretize_reactor(reactor, nfe: int = 10, method: str = 'BACKWARD'):
    """
    Apply discretization to the reactor spatial domain.
    
    Parameters
    ----------
    reactor : SimplifiedRWGSReactorData
        Reactor to discretize
    nfe : int
        Number of finite elements
    method : str
        Discretization method ('BACKWARD' or 'FORWARD')
    """
    discretizer = TransformationFactory('dae.finite_difference')
    discretizer.apply_to(reactor, wrt=reactor.W, nfe=nfe, scheme=method)
    print(f"✓ Reactor discretized with {nfe} finite elements ({method})")


# ==================== SELF-TEST EXAMPLE ====================
if __name__ == "__main__":
    """
    Self-test: Build and solve a minimal RWGS reactor.
    """
    print("\n" + "="*70)
    print("SIMPLIFIED RWGS REACTOR - SELF TEST")
    print("="*70)
    
    from pyomo.environ import ConcreteModel, SolverFactory
    from idaes.core import FlowsheetBlock
    
    # Create model
    print("\n1. Building model...")
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0])
    
    # Define inlet conditions
    inlet_flow = {
        'CO2': 0.7,
        'H2': 0.3,
        'CO': 1e-8,
        'H2O': 1e-8,
    }
    
    # Build reactor - add directly to flowsheet
    print("\n2. Creating reactor...")
    m.fs.reactor = SimplifiedRWGSReactor()
    
    # Initialize
    m.fs.reactor.initialize(
        inlet_flow=inlet_flow,
        temperature=523.15,  # 250°C
        pressure=20e5,  # 20 bar
        W_total=5.0,  # 5 kg catalyst
        k_rwgs=0.05,
        Keq_rwgs=0.8,
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
            print("\n✓ Solution converged!")
            
            # Print results
            t = m.fs.time.first()
            W_inlet = min(m.fs.reactor.W)
            W_outlet = max(m.fs.reactor.W)
            
            print("\n" + "="*70)
            print("RESULTS")
            print("="*70)
            print("\nINLET (W=0):")
            for comp in m.fs.reactor.component_list:
                flow = value(m.fs.reactor.flow_mol_comp[t, W_inlet, comp])
                print(f"  {comp:5s}: {flow:.6f} kmol/s")
            
            print(f"\nOUTLET (W={W_outlet}):")
            for comp in m.fs.reactor.component_list:
                flow = value(m.fs.reactor.flow_mol_comp[t, W_outlet, comp])
                print(f"  {comp:5s}: {flow:.6f} kmol/s")
            
            # Calculate conversions
            CO2_in = value(m.fs.reactor.flow_mol_comp[t, W_inlet, 'CO2'])
            CO2_out = value(m.fs.reactor.flow_mol_comp[t, W_outlet, 'CO2'])
            CO2_conv = (CO2_in - CO2_out) / CO2_in * 100 if CO2_in > 1e-6 else 0
            
            print(f"\nCO2 Conversion: {CO2_conv:.2f}%")
            print("="*70)
            
        else:
            print(f"\n⚠ Solver terminated with: {results.solver.termination_condition}")
    
    except Exception as e:
        print(f"\n✗ Error during solve: {e}")
        import traceback
        traceback.print_exc()
