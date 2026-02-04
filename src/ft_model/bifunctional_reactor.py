"""
Bifunctional Packed-Bed Reactor Unit Model for RWGS + FT + Zeolite Upgrading

This module implements a custom IDAES unit model for a packed-bed reactor
performing three coupled reaction systems:
1. RWGS (Reverse Water-Gas Shift): CO2 + H2 <-> CO + H2O
2. FT (Fischer-Tropsch): CO + H2 -> hydrocarbons
3. Zeolite upgrading: Heavy hydrocarbons -> lighter products

The model uses 1D spatial discretization along the catalyst bed using Pyomo.DAE.
"""
# type: ignore  # Pyomo/IDAES type hints not fully recognized by Pylance

from typing import Dict, List, Optional, Tuple

import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

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

# Simplified FT stoichiometry (CO + 2 H2 consumption)
FT_STOICHIOMETRY = {
    'CH4': {'CO': -1.0, 'H2': -2.0, 'CH4': 1.0},
    'C2H4': {'CO': -2.0, 'H2': -4.0, 'C2H4': 1.0},
    'C2H6': {'CO': -2.0, 'H2': -5.0, 'C2H6': 1.0},
    'C3H6': {'CO': -3.0, 'H2': -6.0, 'C3H6': 1.0},
    'C3H8': {'CO': -3.0, 'H2': -8.0, 'C3H8': 1.0},
    'C5plus_lump': {'CO': -5.0, 'H2': -11.0, 'C5plus_lump': 1.0},
    'wax': {'CO': -20.0, 'H2': -41.0, 'wax': 1.0},
}

# Zeolite cracking stoichiometry
ZEOLITE_STOICHIOMETRY = {
    'distillate': {'wax': -1.0, 'distillate': 0.5, 'LPG': 0.3},
    'naphtha': {'distillate': -1.0, 'naphtha': 0.6, 'LPG': 0.2},
    'aromatics': {'light_olefins': -1.0, 'aromatics': 1.0},
}


@declare_process_block_class('BifunctionalPackedBedReactor')
class BifunctionalPackedBedReactorData(UnitModelBlockData):
    """
    Bifunctional packed-bed reactor combining RWGS, FT, and zeolite upgrading.
    
    Reactions:
    (1) RWGS: CO2 + H2 <-> CO + H2O
    (2) FT: CO + H2 -> CH4, C2H4, C2H6, C3H6, C3H8, C5plus, wax
    (3) Zeolite: wax -> distillate, distillate -> naphtha, olefins -> aromatics
    
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
    
    def build(self):
        """
        Build the bifunctional reactor unit model.
        
        Creates spatial domain, state variables, parameters, and coupled equations.
        """
        super().build()
        
        # Define component list (RWGS feeds + FT products + zeolite products)
        self.component_list = [
            # RWGS components
            'CO2', 'H2', 'CO', 'H2O',
            # FT products
            'CH4', 'C2H4', 'C2H6', 'C3H6', 'C3H8', 'C4_lump', 'C5plus_lump',
            # Zeolite products
            'LPG', 'light_olefins', 'naphtha', 'distillate', 'wax', 
            'aromatics', 'oxygenates', 'coke'
        ]
        
        # Define the spatial domain (catalyst weight, normalized to [0, 1])
        self.W = ContinuousSet(bounds=(0, 1.0))
        
        # ==================== PARAMETERS ====================
        self.W_total = pyo.Param(
            initialize=1.0,
            mutable=True,
            doc='Total catalyst weight [kg]',
        )
        
        # RWGS kinetic parameters
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
        
        # FT kinetic parameters (simplified: one rate constant per product)
        self.k_ft = pyo.Param(
            self.flowsheet().time,
            ['CH4', 'C2H4', 'C3H6', 'C5plus_lump'],
            initialize=0.01,
            mutable=True,
            doc='FT reaction rate constants [kmol/(kg_cat·s·Pa)]',
        )
        
        # Zeolite kinetic parameters
        self.k_zeo = pyo.Param(
            self.flowsheet().time,
            ['wax_cracking', 'distillate_cracking', 'olefin_aromatization'],
            initialize=0.001,
            mutable=True,
            doc='Zeolite cracking rate constants [1/s]',
        )
        
        # ==================== STATE VARIABLES ====================
        # Molar flow rates [kmol/s]
        self.flow_mol_comp = pyo.Var(
            self.flowsheet().time,
            self.W,
            self.component_list,
            initialize=0.1,
            bounds=(1e-8, None),
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
            initialize=1000.0,
            bounds=(0, None),
            doc='Partial pressure of component i',
        )
        
        # Reaction rates [kmol/(kg_cat·s)]
        self.rate_rwgs = pyo.Var(
            self.flowsheet().time,
            self.W,
            initialize=0.001,
            bounds=None,
            doc='RWGS reaction rate',
        )
        
        self.rate_ft = pyo.Var(
            self.flowsheet().time,
            self.W,
            ['CH4', 'C2H4', 'C3H6', 'C5plus_lump'],
            initialize=0.0001,
            bounds=None,
            doc='FT reaction rates for each product',
        )
        
        self.rate_zeo = pyo.Var(
            self.flowsheet().time,
            self.W,
            ['wax_cracking', 'distillate_cracking', 'olefin_aromatization'],
            initialize=0.00001,
            bounds=None,
            doc='Zeolite cracking rates',
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
            y_j = b.flow_mol_comp[t, w, j] / (b.flow_mol[t, w] + 1e-8)
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
                        / (b.Keq_rwgs[t] + 1e-8)
                    )
                )
            )
            return b.rate_rwgs[t, w] == r
        
        # FT reaction rate expressions (simplified, pressure-dependent)
        @self.Constraint(self.flowsheet().time, self.W)
        def rate_ft_ch4_expr(b, t, w):
            """CH4 formation rate from FT."""
            # r_CH4 = k * p_CO * p_H2^2
            r = (
                b.k_ft[t, 'CH4']
                * b.partial_pressure[t, w, 'CO']
                * (b.partial_pressure[t, w, 'H2'] ** 2)
            )
            return b.rate_ft[t, w, 'CH4'] == r
        
        @self.Constraint(self.flowsheet().time, self.W)
        def rate_ft_c2h4_expr(b, t, w):
            """C2H4 formation rate from FT."""
            r = (
                b.k_ft[t, 'C2H4']
                * (b.partial_pressure[t, w, 'CO'] ** 2)
                * (b.partial_pressure[t, w, 'H2'] ** 2)
            )
            return b.rate_ft[t, w, 'C2H4'] == r
        
        @self.Constraint(self.flowsheet().time, self.W)
        def rate_ft_c3h6_expr(b, t, w):
            """C3H6 formation rate from FT."""
            r = (
                b.k_ft[t, 'C3H6']
                * (b.partial_pressure[t, w, 'CO'] ** 3)
                * (b.partial_pressure[t, w, 'H2'] ** 2)
            )
            return b.rate_ft[t, w, 'C3H6'] == r
        
        @self.Constraint(self.flowsheet().time, self.W)
        def rate_ft_c5plus_expr(b, t, w):
            """C5+ and wax formation rate from FT."""
            r = (
                b.k_ft[t, 'C5plus_lump']
                * (b.partial_pressure[t, w, 'CO'] ** 5)
                * (b.partial_pressure[t, w, 'H2'] ** 2)
            )
            return b.rate_ft[t, w, 'C5plus_lump'] == r
        
        # Zeolite cracking rate expressions (first-order in reactant)
        @self.Constraint(self.flowsheet().time, self.W)
        def rate_zeo_wax_cracking(b, t, w):
            """Wax cracking rate."""
            r = b.k_zeo[t, 'wax_cracking'] * b.flow_mol_comp[t, w, 'wax']
            return b.rate_zeo[t, w, 'wax_cracking'] == r
        
        @self.Constraint(self.flowsheet().time, self.W)
        def rate_zeo_distillate_cracking(b, t, w):
            """Distillate cracking rate."""
            r = b.k_zeo[t, 'distillate_cracking'] * b.flow_mol_comp[t, w, 'distillate']
            return b.rate_zeo[t, w, 'distillate_cracking'] == r
        
        @self.Constraint(self.flowsheet().time, self.W)
        def rate_zeo_olefin_arom(b, t, w):
            """Olefin aromatization rate."""
            r = b.k_zeo[t, 'olefin_aromatization'] * b.flow_mol_comp[t, w, 'light_olefins']
            return b.rate_zeo[t, w, 'olefin_aromatization'] == r
        
        # Material balance equations (dF_i/dW = sum_j nu[i,j] * r_j)
        @self.Constraint(
            self.flowsheet().time,
            self.W,
            self.component_list,
        )
        def material_balance(b, t, w, i):
            """
            Material balance for component i.
            
            dF_i/dW = nu_rwgs[i] * r_rwgs 
                    + sum_j nu_ft[i,j] * r_ft[j]
                    + sum_k nu_zeo[i,k] * r_zeo[k]
            """
            # RWGS contribution
            nu_rwgs_i = RWGS_STOICHIOMETRY.get(i, 0.0)
            rwgs_term = nu_rwgs_i * b.rate_rwgs[t, w]
            
            # FT contribution
            ft_term = 0.0
            for ft_product, stoich_dict in FT_STOICHIOMETRY.items():
                nu_ft_ij = stoich_dict.get(i, 0.0)
                ft_term += nu_ft_ij * b.rate_ft[t, w, ft_product]
            
            # Zeolite contribution
            zeo_term = 0.0
            # Wax cracking: wax -> 0.5 distillate + 0.3 LPG
            if i == 'wax':
                zeo_term += -1.0 * b.rate_zeo[t, w, 'wax_cracking']
            elif i == 'distillate':
                zeo_term += 0.5 * b.rate_zeo[t, w, 'wax_cracking']
                zeo_term += -1.0 * b.rate_zeo[t, w, 'distillate_cracking']
            elif i == 'LPG':
                zeo_term += 0.3 * b.rate_zeo[t, w, 'wax_cracking']
                zeo_term += 0.2 * b.rate_zeo[t, w, 'distillate_cracking']
            elif i == 'naphtha':
                zeo_term += 0.6 * b.rate_zeo[t, w, 'distillate_cracking']
            elif i == 'light_olefins':
                zeo_term += -1.0 * b.rate_zeo[t, w, 'olefin_aromatization']
            elif i == 'aromatics':
                zeo_term += 1.0 * b.rate_zeo[t, w, 'olefin_aromatization']
                zeo_term += -0.1 * b.rate_zeo[t, w, 'olefin_aromatization']  # coke from aromatics
            elif i == 'coke':
                zeo_term += 0.1 * b.rate_zeo[t, w, 'olefin_aromatization']
            
            return b.dF_dW[t, w, i] == rwgs_term + ft_term + zeo_term
        
        # Isothermal temperature constraint
        if self.config.isothermal:
            @self.Constraint(self.flowsheet().time, self.W)
            def isothermal_constraint(b, t, w):
                """Temperature is constant along the reactor."""
                return b.temperature[t, w] == b.temperature[t, 0]
        
        # Pressure drop constraint (simplified: constant pressure)
        @self.Constraint(self.flowsheet().time, self.W)
        def pressure_drop_constraint(b, t, w):
            """Assume negligible pressure drop."""
            return b.pressure[t, w] == b.pressure[t, 0]
    
    def initialize(
        self,
        inlet_flow: Dict[str, float],
        inlet_temperature: float,
        inlet_pressure: float,
        W_total_value: float = 1.0,
        k_rwgs_value: float = 0.1,
        Keq_value: float = 1.0,
        k_ft_dict: Optional[Dict[str, float]] = None,
        k_zeo_dict: Optional[Dict[str, float]] = None,
    ):
        """
        Initialize the bifunctional reactor model.
        
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
            RWGS reaction rate constant
        Keq_value : float
            RWGS equilibrium constant
        k_ft_dict : dict, optional
            FT rate constants {'CH4': k, 'C2H4': k, ...}
        k_zeo_dict : dict, optional
            Zeolite rate constants {'wax_cracking': k, ...}
        """
        # Set default kinetic parameters if not provided
        if k_ft_dict is None:
            k_ft_dict = {
                'CH4': 0.01,
                'C2H4': 0.005,
                'C3H6': 0.002,
                'C5plus_lump': 0.0005,
            }
        
        if k_zeo_dict is None:
            k_zeo_dict = {
                'wax_cracking': 0.001,
                'distillate_cracking': 0.0005,
                'olefin_aromatization': 0.0002,
            }
        
        # Set parameter values
        self.W_total.set_value(W_total_value)
        
        time_set = self.flowsheet().time
        for t in time_set:
            self.k_rwgs[t].set_value(k_rwgs_value)
            self.Keq_rwgs[t].set_value(Keq_value)
            
            for ft_product, k_val in k_ft_dict.items():
                self.k_ft[t, ft_product].set_value(k_val)
            
            for zeo_reaction, k_val in k_zeo_dict.items():
                self.k_zeo[t, zeo_reaction].set_value(k_val)
        
        # Initialize state variables
        total_inlet_flow = sum(
            inlet_flow.get(j, 0.0) for j in self.component_list
        )
        
        for t in time_set:
            # Inlet conditions at W=0
            for j in self.component_list:
                self.flow_mol_comp[t, 0, j].set_value(
                    inlet_flow.get(j, 1e-6)
                )
            
            self.temperature[t, 0].set_value(inlet_temperature)
            self.pressure[t, 0].set_value(inlet_pressure)
            self.flow_mol[t, 0].set_value(total_inlet_flow)
            
            # Propagate to interior points with linear profile
            w_points = list(self.W)
            for w in w_points:
                # Linear profile along W
                alpha = w / max(w_points) if max(w_points) > 0 else 0
                for j in self.component_list:
                    inlet_val = inlet_flow.get(j, 1e-6)
                    # Assume some conversion along W
                    conversion = alpha * 0.5 if j in ['CO2', 'H2'] else alpha * 0.3
                    self.flow_mol_comp[t, w, j].set_value(
                        inlet_val * (1 - conversion) if conversion < 1 else 1e-6
                    )
                
                self.temperature[t, w].set_value(inlet_temperature)
                self.pressure[t, w].set_value(inlet_pressure)
    
    def _get_performance_contents(self, time_point=0) -> Dict:
        """
        Report inlet/outlet performance metrics.
        
        Parameters
        ----------
        time_point : int
            Time index to report (default: 0, first time point)
        
        Returns
        -------
        dict
            Dictionary of performance metrics
        """
        t = list(self.flowsheet().time)[time_point]  # type: ignore
        w_inlet = self.W.first()  # type: ignore
        w_outlet = self.W.last()  # type: ignore
        
        # Compile inlet and outlet conditions
        contents = {
            'Inlet': {
                'Temperature [K]': self.temperature[t, w_inlet].value,  # type: ignore
                'Pressure [Pa]': self.pressure[t, w_inlet].value,  # type: ignore
                'Total Molar Flow [kmol/s]': self.flow_mol[t, w_inlet].value,  # type: ignore
                'CO2 Flow [kmol/s]': self.flow_mol_comp[t, w_inlet, 'CO2'].value,  # type: ignore
                'H2 Flow [kmol/s]': self.flow_mol_comp[t, w_inlet, 'H2'].value,  # type: ignore
            },
            'Outlet': {
                'Temperature [K]': self.temperature[t, w_outlet].value,  # type: ignore
                'Pressure [Pa]': self.pressure[t, w_outlet].value,  # type: ignore
                'Total Molar Flow [kmol/s]': self.flow_mol[t, w_outlet].value,  # type: ignore
                'CO Flow [kmol/s]': self.flow_mol_comp[t, w_outlet, 'CO'].value,  # type: ignore
                'CH4 Flow [kmol/s]': self.flow_mol_comp[t, w_outlet, 'CH4'].value,  # type: ignore
                'C5+ Flow [kmol/s]': self.flow_mol_comp[t, w_outlet, 'C5plus_lump'].value,  # type: ignore
                'Gasoline (naphtha) [kmol/s]': self.flow_mol_comp[t, w_outlet, 'naphtha'].value,  # type: ignore
                'Distillate [kmol/s]': self.flow_mol_comp[t, w_outlet, 'distillate'].value,  # type: ignore
            },
            'Kinetics': {
                'RWGS Rate Constant': self.k_rwgs[t].value,  # type: ignore
                'FT CH4 Rate Constant': self.k_ft[t, 'CH4'].value,  # type: ignore
                'Total Catalyst Weight [kg]': self.W_total.value,  # type: ignore
            },
        }
        
        return contents


def build_bifunctional_reactor(
    flowsheet_block,
    property_package,
    inlet_flow: Dict[str, float],
    inlet_temperature: float,
    inlet_pressure: float,
    W_total: float = 1.0,
) -> BifunctionalPackedBedReactorData:
    """
    Helper function to build and initialize a bifunctional reactor.
    
    Parameters
    ----------
    flowsheet_block : Pyomo Block
        Parent flowsheet block
    property_package : object
        Property package for thermodynamic calculations
    inlet_flow : dict
        Inlet molar flows {component: flow_kmol_s}
    inlet_temperature : float
        Inlet temperature [K]
    inlet_pressure : float
        Inlet pressure [Pa]
    W_total : float
        Total catalyst weight [kg]
    
    Returns
    -------
    BifunctionalPackedBedReactorData
        Initialized reactor model
    
    Example
    -------
    >>> reactor = build_bifunctional_reactor(
    ...     flowsheet, props,
    ...     inlet_flow={'CO2': 1.0, 'H2': 2.0, 'CO': 0.1},
    ...     inlet_temperature=523.15,
    ...     inlet_pressure=101325.0,
    ...     W_total=5.0
    ... )
    >>> reactor.initialize(inlet_flow, inlet_temperature, inlet_pressure, W_total)
    """
    # Create reactor block
    reactor = BifunctionalPackedBedReactorData(
        flowsheet_block,
        config={'property_package': property_package, 'isothermal': True}
    )
    
    # Initialize with provided conditions
    reactor.initialize(
        inlet_flow=inlet_flow,
        inlet_temperature=inlet_temperature,
        inlet_pressure=inlet_pressure,
        W_total_value=W_total,
    )
    
    return reactor


if __name__ == "__main__":
    """
    Example usage of bifunctional reactor.
    
    This example shows how to set up and use the bifunctional reactor
    in a flowsheet context.
    """
    print("=" * 70)
    print("Bifunctional Packed-Bed Reactor (RWGS + FT + Zeolite Upgrading)")
    print("=" * 70)
    
    # Example inlet conditions
    inlet_flow = {
        'CO2': 1.0,
        'H2': 2.0,
        'CO': 0.1,
        'H2O': 0.05,
        'CH4': 0.01,
        'C2H4': 0.0,
        'C2H6': 0.0,
        'C3H6': 0.0,
        'C3H8': 0.0,
        'C4_lump': 0.0,
        'C5plus_lump': 0.0,
        'LPG': 0.0,
        'light_olefins': 0.0,
        'naphtha': 0.0,
        'distillate': 0.0,
        'wax': 0.0,
        'aromatics': 0.0,
        'oxygenates': 0.0,
        'coke': 0.0,
    }
    
    inlet_T = 523.15  # 250°C
    inlet_P = 30 * 101325.0  # 30 bar
    W_total = 5.0  # 5 kg catalyst
    
    print("\nInlet Conditions:")
    print(f"  Temperature: {inlet_T:.1f} K ({inlet_T - 273.15:.1f}°C)")
    print(f"  Pressure: {inlet_P/101325:.1f} bar")
    print(f"  Total catalyst weight: {W_total} kg")
    print(f"  CO2 inlet: {inlet_flow['CO2']:.2f} kmol/s")
    print(f"  H2 inlet: {inlet_flow['H2']:.2f} kmol/s")
    
    print("\nReactor Model Features:")
    print("  - 1D spatial discretization along catalyst weight (W)")
    print("  - RWGS reaction: CO2 + H2 <-> CO + H2O")
    print("  - FT synthesis: CO + H2 -> CH4, C2H4, C5+, wax")
    print("  - Zeolite upgrading: wax -> distillate -> naphtha, olefins -> aromatics")
    print("  - Coupled material balance equations with DerivativeVar")
    print("  - Partial pressure calculations from mole fractions")
    print("  - Isothermal operation constraint")
    
    print("\nComponent List:")
    reactor_data = BifunctionalPackedBedReactorData()
    reactor_data.component_list = [
        'CO2', 'H2', 'CO', 'H2O',
        'CH4', 'C2H4', 'C2H6', 'C3H6', 'C3H8', 'C4_lump', 'C5plus_lump',
        'LPG', 'light_olefins', 'naphtha', 'distillate', 'wax',
        'aromatics', 'oxygenates', 'coke'
    ]
    for i, comp in enumerate(reactor_data.component_list, 1):
        print(f"  {i:2d}. {comp}")
    
    print("\n" + "=" * 70)
    print("To use in a flowsheet:")
    print("=" * 70)
    print("""
    from src.ft_model.bifunctional_reactor import build_bifunctional_reactor
    
    # Build flowsheet
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0])
    
    # Add property package
    m.fs.properties = PropertyPackageClass()
    
    # Build reactor
    m.fs.reactor = build_bifunctional_reactor(
        m.fs, m.fs.properties,
        inlet_flow=inlet_flow,
        inlet_temperature=inlet_T,
        inlet_pressure=inlet_P,
        W_total=W_total
    )
    
    # Connect feed streams and solve
    # ...
    """)
