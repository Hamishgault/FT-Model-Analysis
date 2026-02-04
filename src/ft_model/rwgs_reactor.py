"""
Packed-Bed Plug Flow Reactor Unit Model for Reverse Water-Gas Shift (RWGS)

Reaction: CO2 + H2 <-> CO + H2O

This module implements a custom IDAES unit model for a packed-bed reactor
with 1D spatial discretization along the catalyst bed using Pyomo.DAE.
"""
# type: ignore  # Pyomo/IDAES type hints not fully recognized by Pylance

from typing import Dict, List, Optional

import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

from idaes.core import (
    UnitModelBlockData,
    declare_process_block_class,
)


# Stoichiometric coefficients for RWGS reaction: CO2 + H2 <-> CO + H2O
RWGS_STOICHIOMETRY = {
    'CO2': -1.0,   # consumed
    'H2': -1.0,    # consumed
    'CO': 1.0,     # produced
    'H2O': 1.0,    # produced
}


@declare_process_block_class('PackedBedRWGSReactor')
class PackedBedRWGSReactorData(UnitModelBlockData):
    """
    Packed-bed plug flow reactor for reverse water-gas shift (RWGS) reaction.
    
    Reaction: CO2 + H2 <-> CO + H2O
    
    This model uses 1D spatial discretization along the catalyst weight (W) dimension.
    
    Attributes
    ----------
    inlet : Port
        Inlet stream to the reactor
    outlet : Port
        Outlet stream from the reactor
    W : ContinuousSet
        Catalyst weight domain [0, W_total]
    flow_mol_comp : Var
        Molar flow rate of each component [kmol/s] indexed by (time, W, component)
    temperature : Var
        Temperature along the reactor [K] indexed by (time, W)
    pressure : Var
        Pressure along the reactor [Pa] indexed by (time, W)
    """
    
    CONFIG = UnitModelBlockData.CONFIG()
    # Configuration options are inherited from UnitModelBlockData
    # Subclasses can override specific settings as needed
    
    
    def build(self):
        """
        Build the unit model.
        
        Creates spatial domain, state variables, parameters, and equations.
        """
        super().build()
        
        # Define component list for RWGS reaction
        self.component_list = ['CO2', 'H2', 'CO', 'H2O']
        
        # Define the spatial domain (catalyst weight in kg)
        self.W = ContinuousSet(bounds=(0, 1.0))  # Normalized to [0, 1], actual weight is W_total * W
        
        # ==================== PARAMETERS ====================
        self.W_total = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Total catalyst weight [kg]',
        )
        
        self.k_rwgs = pyo.Param(
            self.flowsheet().time,
            initialize=0.1,
            mutable=True,
            doc='RWGS reaction rate constant [kmol/(kg_cat·s·Pa^2)]',
        )
        
        self.Keq_rwgs = pyo.Param(
            self.flowsheet().time,
            initialize=1.0,
            mutable=True,
            doc='RWGS equilibrium constant (dimensionless)',
        )
        
        # ==================== STATE VARIABLES ====================
        # Molar flow rates [kmol/s]
        self.flow_mol_comp = pyo.Var(
            self.flowsheet().time,
            self.W,
            self.component_list,
            initialize=1.0,
            bounds=(1e-6, None),
            doc='Molar flow rate of component i',
        )
        
        # Temperature [K]
        self.temperature = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=500.0,
            bounds=(200.0, 1000.0),
            doc='Temperature along the reactor',
        )
        
        # Pressure [Pa]
        self.pressure = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=101325.0,
            bounds=(50000, 3000000),
            doc='Pressure along the reactor',
        )
        
        # Total molar flow [kmol/s]
        self.flow_mol = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=1.0,
            bounds=(1e-6, None),
            doc='Total molar flow rate',
        )
        
        # ==================== DERIVED VARIABLES ====================
        # Partial pressures [Pa]
        self.partial_pressure = pyo.Var(
            self.flowsheet().time,
            self.W,
            self.component_list,
            initialize=101325.0,
            bounds=(0, None),
            doc='Partial pressure of component i',
        )
        
        # Reaction rate [kmol/(kg_cat·s)]
        self.rate_rwgs = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.01,
            bounds=None,
            doc='RWGS reaction rate',
        )
        
        # ==================== DIFFERENTIAL EQUATIONS ====================
        # Material balance derivatives
        self.dF_dW = DerivativeVar(
            self.flow_mol_comp,
            wrt=self.W,
            initialize=0.0,
            doc='Material balance derivatives',
        )
        
        # ==================== CONSTRAINTS ====================
        # Total molar flow definition
        @self.Constraint(self.flowsheet().time, self.W)
        def total_flow_balance(b, t, w):
            """Total flow is sum of component flows."""
            return b.flow_mol[t, w] == sum(
                b.flow_mol_comp[t, w, j] for j in b.component_list
            )
        
        # Partial pressure calculation
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            self.component_list,
        )
        def partial_pressure_calc(b, t, w, j):
            """Calculate partial pressure from mole fraction and total pressure."""
            y_j = b.flow_mol_comp[t, w, j] / b.flow_mol[t, w]
            return b.partial_pressure[t, w, j] == y_j * b.pressure[t, w]
        
        # RWGS reaction rate expression
        @self.Constraint(self.flowsheet().time, self.W)
        def rate_rwgs_expr(b, t, w):
            """
            RWGS rate expression.
            
            r = k * (p_CO2 * p_H2 - p_CO * p_H2O / Keq)
            """
            r = (
                b.k_rwgs[t]
                * (
                    b.partial_pressure[t, w, 'CO2']
                    * b.partial_pressure[t, w, 'H2']
                    - (
                        b.partial_pressure[t, w, 'CO']
                        * b.partial_pressure[t, w, 'H2O']
                        / (b.Keq_rwgs[t] + 1e-8)  # small offset to avoid division by zero
                    )
                )
            )
            return b.rate_rwgs[t, w] == r
        
        # Material balance equations (dF_i/dW = nu_i * r)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            self.component_list,
        )
        def material_balance(b, t, w, j):
            """
            Material balance for component j.
            
            dF_j/dW = nu_j * r_rwgs
            """
            nu_j = RWGS_STOICHIOMETRY.get(j, 0.0)  # 0 for non-reactive species
            return b.dF_dW[t, w, j] == nu_j * b.rate_rwgs[t, w]
        
        # Isothermal temperature constraint (if isothermal flag is True)
        if self.config.isothermal:
            @self.Constraint(self.flowsheet().time, self.W)
            def isothermal_constraint(b, t, w):
                """Temperature is constant along the reactor."""
                return b.temperature[t, w] == b.temperature[t, 0]
        
        # Pressure drop (optional: simple linear drop assumption)
        @self.Constraint(self.flowsheet().time, self.W)
        def pressure_drop_constraint(b, t, w):
            """
            Simple linear pressure drop along the reactor.
            P(W) = P_inlet - dP/dW * W * W_total
            where dP/dW is a small fraction of inlet pressure.
            """
            # For now, assume negligible pressure drop
            return b.pressure[t, w] == b.pressure[t, 0]
    
    def initialize(
        self,
        inlet_flow: Dict[str, float],
        inlet_temperature: float,
        inlet_pressure: float,
        W_total_value: float = 1.0,
        k_rwgs_value: float = 0.1,
        Keq_value: float = 1.0,
    ):
        """
        Initialize the reactor model.
        
        Parameters
        ----------
        inlet_flow : dict
            Inlet molar flows {component: flow_kmol_s}
        inlet_temperature : float
            Inlet temperature [K]
        inlet_pressure : float
            Inlet pressure [Pa]
        W_total_value : float
            Total catalyst weight [kg]
        k_rwgs_value : float
            Reaction rate constant
        Keq_value : float
            Equilibrium constant
        """
        props = self.config.property_package
        
        # Set parameter values
        self.W_total.set_value(W_total_value)
        
        time_set = self.flowsheet().time
        for t in time_set:
            self.k_rwgs[t].set_value(k_rwgs_value)
            self.Keq_rwgs[t].set_value(Keq_value)
        
        # Initialize state variables at W=0 (inlet)
        total_inlet_flow = sum(inlet_flow.values())
        
        for t in time_set:
            # Set inlet conditions at W=0
            for j in self.component_list:
                self.flow_mol_comp[t, 0, j].set_value(
                    inlet_flow.get(j, 1e-6)
                )
            
            self.flow_mol[t, 0].set_value(total_inlet_flow)
            self.temperature[t, 0].set_value(inlet_temperature)
            self.pressure[t, 0].set_value(inlet_pressure)
            
            # Propagate initial conditions along W with simple conversion
            # Assume linear conversion profile as initial guess
            for w in self.W:
                if w > 0:
                    # Simple linear conversion guess
                    conversion_frac = min(w, 0.5)  # Assume at most 50% conversion
                    
                    for j in self.component_list:
                        inlet_flow_j = inlet_flow.get(j, 1e-6)
                        nu_j = RWGS_STOICHIOMETRY.get(j, 0.0)
                        
                        # Update flow based on stoichiometry
                        if 'CO2' in inlet_flow and inlet_flow['CO2'] > 0:
                            extent = conversion_frac * inlet_flow['CO2']
                            self.flow_mol_comp[t, w, j].set_value(
                                inlet_flow_j + nu_j * extent
                            )
                        else:
                            self.flow_mol_comp[t, w, j].set_value(inlet_flow_j)
                    
                    # Update total flow
                    self.flow_mol[t, w].set_value(
                        sum(
                            self.flow_mol_comp[t, w, j].value
                            for j in self.component_list
                        )
                    )
                    
                    self.temperature[t, w].set_value(inlet_temperature)
                    self.pressure[t, w].set_value(inlet_pressure)
    
    def _get_performance_contents(self, time_point=0):
        """
        Return performance contents for reporting.
        
        Parameters
        ----------
        time_point : int
            Time index for reporting
        
        Returns
        -------
        dict
            Dictionary of performance metrics
        """
        t = list(self.flowsheet().time)[time_point]
        w_inlet = self.W.first()
        w_outlet = self.W.last()
        
        perf_dict = {
            'Inlet Conditions': {
                'Temperature [K]': self.temperature[t, w_inlet].value,
                'Pressure [Pa]': self.pressure[t, w_inlet].value,
                'Total Molar Flow [kmol/s]': self.flow_mol[t, w_inlet].value,
            },
            'Outlet Conditions': {
                'Temperature [K]': self.temperature[t, w_outlet].value,
                'Pressure [Pa]': self.pressure[t, w_outlet].value,
                'Total Molar Flow [kmol/s]': self.flow_mol[t, w_outlet].value,
            },
            'Reaction Parameters': {
                'Rate Constant k_rwgs': self.k_rwgs[t].value,
                'Equilibrium Constant Keq': self.Keq_rwgs[t].value,
                'Total Catalyst Weight [kg]': self.W_total.value,
            },
        }
        
        return perf_dict


class PackedBedRWGSReactor(PackedBedRWGSReactorData):
    """
    Packed-bed plug flow reactor for RWGS reaction.
    
    This is the main class that users interact with.
    """
    pass
