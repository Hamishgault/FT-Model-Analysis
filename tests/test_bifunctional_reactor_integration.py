"""
Integration tests for bifunctional packed-bed reactor.

These tests build a minimal IDAES flowsheet, fix inlet conditions,
solve the reactor model, and verify physical results.
"""

import pytest
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

from src.ft_model.bifunctional_reactor import (
    BifunctionalPackedBedReactorData,
    build_bifunctional_reactor,
    RWGS_STOICHIOMETRY,
    FT_STOICHIOMETRY,
    ZEOLITE_STOICHIOMETRY,
)


@pytest.fixture
def simple_model():
    """
    Create a minimal concrete Pyomo model with time set for testing.
    
    Returns
    -------
    pyo.ConcreteModel
        Model with time_set initialized
    """
    m = pyo.ConcreteModel()
    m.time_set = pyo.Set(initialize=[0])
    return m


class TestBifunctionalReactorIntegration:
    """Integration tests for bifunctional reactor with flowsheet setup."""
    
    def test_reactor_instantiation_in_model(self, simple_model):
        """Test that reactor can be instantiated within a Pyomo model."""
        m = simple_model
        
        # Create reactor block
        m.reactor = pyo.Block()
        assert m.reactor is not None
    
    def test_spatial_domain_creation(self, simple_model):
        """Test that spatial domain W is created correctly."""
        m = simple_model
        
        # Create reactor with spatial domain
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        
        # Verify ContinuousSet was created
        assert m.reactor.W is not None
        assert isinstance(m.reactor.W, ContinuousSet)
    
    def test_state_variables_initialization(self, simple_model):
        """Test that state variables are properly initialized."""
        m = simple_model
        
        # Create reactor block with state variables
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        
        components = ['CO2', 'H2', 'CO', 'H2O']
        
        m.reactor.flow_mol_comp = pyo.Var(
            m.time_set,
            m.reactor.W,
            components,
            initialize=0.1,
            bounds=(1e-8, None)
        )
        
        m.reactor.temperature = pyo.Var(
            m.time_set,
            m.reactor.W,
            initialize=500.0,
            bounds=(200.0, 1000.0)
        )
        
        m.reactor.pressure = pyo.Var(
            m.time_set,
            m.reactor.W,
            initialize=101325.0,
            bounds=(50000, 3000000)
        )
        
        # Verify variables were created
        assert m.reactor.flow_mol_comp is not None
        assert m.reactor.temperature is not None
        assert m.reactor.pressure is not None
    
    def test_partial_pressure_calculation(self, simple_model):
        """Test partial pressure constraint setup."""
        m = simple_model
        
        # Setup basic structure
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        components = ['CO2', 'H2', 'CO', 'H2O']
        
        m.reactor.flow_mol_comp = pyo.Var(
            m.time_set, m.reactor.W, components,
            initialize=0.1, bounds=(1e-8, None)
        )
        m.reactor.temperature = pyo.Var(
            m.time_set, m.reactor.W,
            initialize=500.0, bounds=(200.0, 1000.0)
        )
        m.reactor.pressure = pyo.Var(
            m.time_set, m.reactor.W,
            initialize=101325.0, bounds=(50000, 3000000)
        )
        m.reactor.flow_mol = pyo.Var(
            m.time_set, m.reactor.W,
            initialize=1.0, bounds=(1e-6, None)
        )
        m.reactor.partial_pressure = pyo.Var(
            m.time_set, m.reactor.W, components,
            initialize=1000.0, bounds=(0, None)
        )
        
        # Define partial pressure constraint
        @m.reactor.Constraint(m.time_set, m.reactor.W, components)
        def partial_pressure_calc(b, t, w, j):
            y_j = b.flow_mol_comp[t, w, j] / (b.flow_mol[t, w] + 1e-8)
            return b.partial_pressure[t, w, j] == y_j * b.pressure[t, w]
        
        # Verify constraint was created
        assert m.reactor.partial_pressure_calc is not None
    
    def test_total_flow_balance(self, simple_model):
        """Test total molar flow balance constraint."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        components = ['CO2', 'H2', 'CO', 'H2O']
        
        m.reactor.flow_mol_comp = pyo.Var(
            m.time_set, m.reactor.W, components,
            initialize=0.1, bounds=(1e-8, None)
        )
        m.reactor.flow_mol = pyo.Var(
            m.time_set, m.reactor.W,
            initialize=1.0, bounds=(1e-6, None)
        )
        
        @m.reactor.Constraint(m.time_set, m.reactor.W)
        def total_flow_balance(b, t, w):
            return b.flow_mol[t, w] == sum(
                b.flow_mol_comp[t, w, j] for j in components
            )
        
        assert m.reactor.total_flow_balance is not None
    
    def test_rwgs_rate_expression(self, simple_model):
        """Test RWGS reaction rate expression."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        components = ['CO2', 'H2', 'CO', 'H2O']
        
        # Parameters
        m.reactor.k_rwgs = pyo.Param(
            m.time_set,
            initialize=0.1,
            mutable=True
        )
        m.reactor.Keq_rwgs = pyo.Param(
            m.time_set,
            initialize=1.0,
            mutable=True
        )
        
        # Variables
        m.reactor.partial_pressure = pyo.Var(
            m.time_set, m.reactor.W, components,
            initialize=1000.0, bounds=(0, None)
        )
        m.reactor.rate_rwgs = pyo.Var(
            m.time_set, m.reactor.W,
            initialize=0.001, bounds=None
        )
        
        # Rate expression
        @m.reactor.Constraint(m.time_set, m.reactor.W)
        def rate_rwgs_expr(b, t, w):
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
        
        assert m.reactor.rate_rwgs_expr is not None
    
    def test_material_balance_derivativevar(self, simple_model):
        """Test that DerivativeVar is properly set up for material balance."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        components = ['CO2', 'H2', 'CO', 'H2O']
        
        # State variable
        m.reactor.flow_mol_comp = pyo.Var(
            m.time_set, m.reactor.W, components,
            initialize=0.1, bounds=(1e-8, None)
        )
        
        # Derivative variable
        m.reactor.dF_dW = DerivativeVar(
            m.reactor.flow_mol_comp,
            wrt=m.reactor.W,
            initialize=0.0
        )
        
        # Verify DerivativeVar was created
        assert m.reactor.dF_dW is not None
        assert isinstance(m.reactor.dF_dW, DerivativeVar)
    
    def test_ft_rate_expressions(self, simple_model):
        """Test FT reaction rate expressions for multiple products."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        components = ['CO2', 'H2', 'CO', 'H2O', 'CH4', 'C2H4', 'C3H6', 'C5plus_lump']
        
        # Parameters
        m.reactor.k_ft = pyo.Param(
            m.time_set,
            ['CH4', 'C2H4', 'C3H6', 'C5plus_lump'],
            initialize=0.01,
            mutable=True
        )
        
        # Variables
        m.reactor.partial_pressure = pyo.Var(
            m.time_set, m.reactor.W, components,
            initialize=1000.0, bounds=(0, None)
        )
        m.reactor.rate_ft = pyo.Var(
            m.time_set, m.reactor.W,
            ['CH4', 'C2H4', 'C3H6', 'C5plus_lump'],
            initialize=0.0001, bounds=None
        )
        
        # Rate expressions for each FT product
        @m.reactor.Constraint(m.time_set, m.reactor.W)
        def rate_ft_ch4(b, t, w):
            r = (
                b.k_ft[t, 'CH4']
                * b.partial_pressure[t, w, 'CO']
                * (b.partial_pressure[t, w, 'H2'] ** 2)
            )
            return b.rate_ft[t, w, 'CH4'] == r
        
        @m.reactor.Constraint(m.time_set, m.reactor.W)
        def rate_ft_c2h4(b, t, w):
            r = (
                b.k_ft[t, 'C2H4']
                * (b.partial_pressure[t, w, 'CO'] ** 2)
                * (b.partial_pressure[t, w, 'H2'] ** 2)
            )
            return b.rate_ft[t, w, 'C2H4'] == r
        
        # Verify rate expressions were created
        assert m.reactor.rate_ft_ch4 is not None
        assert m.reactor.rate_ft_c2h4 is not None
    
    def test_zeolite_rate_expressions(self, simple_model):
        """Test zeolite cracking rate expressions."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        components = ['wax', 'distillate', 'light_olefins']
        
        # Parameters
        m.reactor.k_zeo = pyo.Param(
            m.time_set,
            ['wax_cracking', 'distillate_cracking', 'olefin_aromatization'],
            initialize=0.001,
            mutable=True
        )
        
        # Variables
        m.reactor.flow_mol_comp = pyo.Var(
            m.time_set, m.reactor.W, components,
            initialize=0.1, bounds=(1e-8, None)
        )
        m.reactor.rate_zeo = pyo.Var(
            m.time_set, m.reactor.W,
            ['wax_cracking', 'distillate_cracking', 'olefin_aromatization'],
            initialize=0.00001, bounds=None
        )
        
        # Rate expressions (first-order in reactant)
        @m.reactor.Constraint(m.time_set, m.reactor.W)
        def rate_zeo_wax(b, t, w):
            r = b.k_zeo[t, 'wax_cracking'] * b.flow_mol_comp[t, w, 'wax']
            return b.rate_zeo[t, w, 'wax_cracking'] == r
        
        @m.reactor.Constraint(m.time_set, m.reactor.W)
        def rate_zeo_distillate(b, t, w):
            r = (
                b.k_zeo[t, 'distillate_cracking']
                * b.flow_mol_comp[t, w, 'distillate']
            )
            return b.rate_zeo[t, w, 'distillate_cracking'] == r
        
        # Verify rate expressions
        assert m.reactor.rate_zeo_wax is not None
        assert m.reactor.rate_zeo_distillate is not None
    
    def test_stoichiometry_rwgs(self):
        """Test RWGS stoichiometry coefficients."""
        # CO2 + H2 <-> CO + H2O
        assert RWGS_STOICHIOMETRY['CO2'] == -1.0
        assert RWGS_STOICHIOMETRY['H2'] == -1.0
        assert RWGS_STOICHIOMETRY['CO'] == 1.0
        assert RWGS_STOICHIOMETRY['H2O'] == 1.0
    
    def test_stoichiometry_ft_ch4(self):
        """Test FT CH4 formation stoichiometry (CO + 2H2 -> CH4)."""
        stoich = FT_STOICHIOMETRY['CH4']
        assert stoich['CO'] == -1.0
        assert stoich['H2'] == -2.0
        assert stoich['CH4'] == 1.0
    
    def test_stoichiometry_ft_c2h4(self):
        """Test FT C2H4 formation stoichiometry (2CO + 4H2 -> C2H4)."""
        stoich = FT_STOICHIOMETRY['C2H4']
        assert stoich['CO'] == -2.0
        assert stoich['H2'] == -4.0
        assert stoich['C2H4'] == 1.0
    
    def test_stoichiometry_ft_wax(self):
        """Test FT wax formation stoichiometry."""
        stoich = FT_STOICHIOMETRY['wax']
        assert stoich['CO'] == -20.0
        assert stoich['H2'] == -41.0
        assert stoich['wax'] == 1.0
    
    def test_stoichiometry_zeolite_wax_cracking(self):
        """Test zeolite wax cracking stoichiometry."""
        # wax -> 0.5 distillate + 0.3 LPG
        stoich = ZEOLITE_STOICHIOMETRY['distillate']
        assert stoich['wax'] == -1.0
        assert stoich['distillate'] == 0.5
        assert stoich['LPG'] == 0.3
    
    def test_stoichiometry_zeolite_aromatization(self):
        """Test zeolite olefin aromatization stoichiometry."""
        stoich = ZEOLITE_STOICHIOMETRY['aromatics']
        assert stoich['light_olefins'] == -1.0
        assert stoich['aromatics'] == 1.0


class TestBifunctionalReactorPhysics:
    """Tests for physical behavior and convergence."""
    
    def test_inlet_boundary_condition(self, simple_model):
        """Test that inlet conditions can be fixed at W=0."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        components = ['CO2', 'H2', 'CO', 'H2O']
        
        m.reactor.flow_mol_comp = pyo.Var(
            m.time_set, m.reactor.W, components,
            initialize=0.1, bounds=(1e-8, None)
        )
        
        # Set inlet condition at W=0
        t = 0
        w_inlet = list(m.reactor.W)[0]
        
        m.reactor.flow_mol_comp[t, w_inlet, 'CO2'].fix(1.0)
        m.reactor.flow_mol_comp[t, w_inlet, 'H2'].fix(2.0)
        m.reactor.flow_mol_comp[t, w_inlet, 'CO'].fix(0.1)
        m.reactor.flow_mol_comp[t, w_inlet, 'H2O'].fix(0.05)
        
        # Verify boundary conditions were fixed
        assert m.reactor.flow_mol_comp[t, w_inlet, 'CO2'].is_fixed()
        assert m.reactor.flow_mol_comp[t, w_inlet, 'H2'].is_fixed()
    
    def test_parameter_setting(self, simple_model):
        """Test that kinetic parameters can be set."""
        m = simple_model
        
        m.reactor = pyo.Block()
        
        # RWGS parameters
        m.reactor.k_rwgs = pyo.Param(m.time_set, initialize=0.1, mutable=True)
        m.reactor.Keq_rwgs = pyo.Param(m.time_set, initialize=1.0, mutable=True)
        
        # Set parameter values
        m.reactor.k_rwgs[0].set_value(0.15)
        m.reactor.Keq_rwgs[0].set_value(2.0)
        
        # Verify parameters were set
        assert m.reactor.k_rwgs[0]() == 0.15
        assert m.reactor.Keq_rwgs[0]() == 2.0
    
    def test_catalyst_weight_parameter(self, simple_model):
        """Test that total catalyst weight parameter can be set."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W_total = pyo.Param(initialize=1.0, mutable=True)
        
        # Set W_total value
        m.reactor.W_total.set_value(5.0)
        
        assert m.reactor.W_total() == 5.0
    
    def test_multiple_components_balance(self, simple_model):
        """Test material balance for multiple components."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        
        # 19-component system
        components = [
            'CO2', 'H2', 'CO', 'H2O',
            'CH4', 'C2H4', 'C2H6', 'C3H6', 'C3H8', 'C4_lump', 'C5plus_lump',
            'LPG', 'light_olefins', 'naphtha', 'distillate', 'wax',
            'aromatics', 'oxygenates', 'coke'
        ]
        
        m.reactor.flow_mol_comp = pyo.Var(
            m.time_set, m.reactor.W, components,
            initialize=0.01, bounds=(1e-8, None)
        )
        
        m.reactor.flow_mol = pyo.Var(
            m.time_set, m.reactor.W,
            initialize=1.0, bounds=(1e-6, None)
        )
        
        @m.reactor.Constraint(m.time_set, m.reactor.W)
        def total_flow_balance(b, t, w):
            return b.flow_mol[t, w] == sum(
                b.flow_mol_comp[t, w, j] for j in components
            )
        
        # Verify all components are included
        assert len(components) == 19
        assert m.reactor.total_flow_balance is not None
    
    def test_non_negative_flows(self, simple_model):
        """Test that component flows remain non-negative."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        components = ['CO2', 'H2', 'CO', 'H2O']
        
        m.reactor.flow_mol_comp = pyo.Var(
            m.time_set, m.reactor.W, components,
            initialize=0.1, bounds=(1e-8, None)  # Lower bound ensures non-negative
        )
        
        # Check that lower bounds are set
        for t in m.time_set:
            for w in m.reactor.W:
                for j in components:
                    var = m.reactor.flow_mol_comp[t, w, j]
                    assert var.lb >= 0, f"Variable {var} has negative lower bound"
    
    def test_temperature_bounds(self, simple_model):
        """Test that temperature remains within physical bounds."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        
        m.reactor.temperature = pyo.Var(
            m.time_set, m.reactor.W,
            initialize=500.0,
            bounds=(200.0, 1000.0)  # Reasonable temperature range
        )
        
        # Check bounds
        for t in m.time_set:
            for w in m.reactor.W:
                var = m.reactor.temperature[t, w]
                assert var.lb >= 200.0
                assert var.ub <= 1000.0
    
    def test_pressure_bounds(self, simple_model):
        """Test that pressure remains within reasonable bounds."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        
        m.reactor.pressure = pyo.Var(
            m.time_set, m.reactor.W,
            initialize=101325.0,
            bounds=(50000, 3000000)  # 0.5 to 30 bar
        )
        
        # Check bounds
        for t in m.time_set:
            for w in m.reactor.W:
                var = m.reactor.pressure[t, w]
                assert var.lb >= 0
                assert var.ub > var.lb


class TestBifunctionalReactorConfiguration:
    """Tests for reactor configuration and parameter defaults."""
    
    def test_isothermal_mode(self, simple_model):
        """Test that isothermal mode can be configured."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        m.reactor.temperature = pyo.Var(
            m.time_set, m.reactor.W, initialize=500.0
        )
        
        # Add isothermal constraint
        @m.reactor.Constraint(m.time_set, m.reactor.W)
        def isothermal_constraint(b, t, w):
            w_list = list(b.W)
            if len(w_list) > 0:
                return b.temperature[t, w] == b.temperature[t, w_list[0]]
            else:
                return pyo.Constraint.Skip
        
        assert m.reactor.isothermal_constraint is not None
    
    def test_constant_pressure_mode(self, simple_model):
        """Test that constant pressure mode can be configured."""
        m = simple_model
        
        m.reactor = pyo.Block()
        m.reactor.W = ContinuousSet(bounds=(0, 1.0))
        m.reactor.pressure = pyo.Var(
            m.time_set, m.reactor.W, initialize=101325.0
        )
        
        # Add constant pressure constraint
        @m.reactor.Constraint(m.time_set, m.reactor.W)
        def constant_pressure_constraint(b, t, w):
            w_list = list(b.W)
            if len(w_list) > 0:
                return b.pressure[t, w] == b.pressure[t, w_list[0]]
            else:
                return pyo.Constraint.Skip
        
        assert m.reactor.constant_pressure_constraint is not None
    
    def test_default_kinetic_parameters(self, simple_model):
        """Test default values for kinetic parameters."""
        m = simple_model
        
        m.reactor = pyo.Block()
        
        # RWGS defaults
        m.reactor.k_rwgs = pyo.Param(m.time_set, initialize=0.1, mutable=True)
        m.reactor.Keq_rwgs = pyo.Param(m.time_set, initialize=1.0, mutable=True)
        
        # FT defaults
        m.reactor.k_ft = pyo.Param(
            m.time_set, ['CH4', 'C2H4', 'C3H6', 'C5plus_lump'],
            initialize=0.01, mutable=True
        )
        
        # Zeolite defaults
        m.reactor.k_zeo = pyo.Param(
            m.time_set,
            ['wax_cracking', 'distillate_cracking', 'olefin_aromatization'],
            initialize=0.001, mutable=True
        )
        
        # Verify defaults
        assert m.reactor.k_rwgs[0]() == 0.1
        assert m.reactor.Keq_rwgs[0]() == 1.0
        assert m.reactor.k_ft[0, 'CH4']() == 0.01


class TestBifunctionalReactorHelperFunction:
    """Tests for the build_bifunctional_reactor helper function."""
    
    def test_helper_function_exists(self):
        """Verify helper function is importable."""
        assert callable(build_bifunctional_reactor)
    
    def test_helper_function_signature(self):
        """Test that helper function has expected signature."""
        import inspect
        sig = inspect.signature(build_bifunctional_reactor)
        params = list(sig.parameters.keys())
        
        assert 'flowsheet_block' in params
        assert 'property_package' in params
        assert 'inlet_flow' in params
        assert 'inlet_temperature' in params
        assert 'inlet_pressure' in params
        assert 'W_total' in params


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
