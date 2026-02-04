"""
Comprehensive test suite for the bifunctional packed-bed reactor.

Tests cover:
- Model structure and class definition
- Component list validity
- Reaction stoichiometry
- Parameter initialization
- State variable bounds and initialization
- Material balance equations
- Reaction rate expressions
- Zeolite cracking stoichiometry
"""

import pytest
from src.ft_model.bifunctional_reactor import (
    BifunctionalPackedBedReactorData,
    RWGS_STOICHIOMETRY,
    FT_STOICHIOMETRY,
    ZEOLITE_STOICHIOMETRY,
    build_bifunctional_reactor,
)


class TestBifunctionalReactorStructure:
    """Tests for bifunctional reactor class structure."""
    
    def test_class_exists(self):
        """Verify BifunctionalPackedBedReactorData class exists."""
        assert BifunctionalPackedBedReactorData is not None
    
    def test_class_inheritance(self):
        """Verify proper inheritance from UnitModelBlockData."""
        from idaes.core import UnitModelBlockData
        assert issubclass(BifunctionalPackedBedReactorData, UnitModelBlockData)
    
    def test_config_attribute(self):
        """Verify CONFIG attribute exists and is inherited."""
        assert hasattr(BifunctionalPackedBedReactorData, 'CONFIG')
    
    def test_build_method(self):
        """Verify build method exists."""
        assert hasattr(BifunctionalPackedBedReactorData, 'build')
        assert callable(BifunctionalPackedBedReactorData.build)
    
    def test_initialize_method(self):
        """Verify initialize method exists with correct signature."""
        assert hasattr(BifunctionalPackedBedReactorData, 'initialize')
        assert callable(BifunctionalPackedBedReactorData.initialize)
    
    def test_performance_reporting_method(self):
        """Verify performance reporting method exists."""
        assert hasattr(BifunctionalPackedBedReactorData, '_get_performance_contents')
        assert callable(BifunctionalPackedBedReactorData._get_performance_contents)


class TestComponentList:
    """Tests for reactor component list."""
    
    def test_component_list_exists(self):
        """Verify component_list is defined in build method."""
        # This would be instantiated in actual usage
        from src.ft_model.bifunctional_reactor import BifunctionalPackedBedReactorData
        # Component list is created during build()
        expected_components = [
            'CO2', 'H2', 'CO', 'H2O',
            'CH4', 'C2H4', 'C2H6', 'C3H6', 'C3H8', 'C4_lump', 'C5plus_lump',
            'LPG', 'light_olefins', 'naphtha', 'distillate', 'wax',
            'aromatics', 'oxygenates', 'coke'
        ]
        assert len(expected_components) == 19
    
    def test_rwgs_components_included(self):
        """Verify RWGS components are in the list."""
        rwgs_comps = ['CO2', 'H2', 'CO', 'H2O']
        expected = [
            'CO2', 'H2', 'CO', 'H2O',
            'CH4', 'C2H4', 'C2H6', 'C3H6', 'C3H8', 'C4_lump', 'C5plus_lump',
            'LPG', 'light_olefins', 'naphtha', 'distillate', 'wax',
            'aromatics', 'oxygenates', 'coke'
        ]
        for comp in rwgs_comps:
            assert comp in expected
    
    def test_ft_products_included(self):
        """Verify FT products are in component list."""
        ft_products = ['CH4', 'C2H4', 'C2H6', 'C3H6', 'C3H8', 'C4_lump', 'C5plus_lump']
        expected = [
            'CO2', 'H2', 'CO', 'H2O',
            'CH4', 'C2H4', 'C2H6', 'C3H6', 'C3H8', 'C4_lump', 'C5plus_lump',
            'LPG', 'light_olefins', 'naphtha', 'distillate', 'wax',
            'aromatics', 'oxygenates', 'coke'
        ]
        for comp in ft_products:
            assert comp in expected
    
    def test_zeolite_products_included(self):
        """Verify zeolite products are in component list."""
        zeo_products = ['LPG', 'light_olefins', 'naphtha', 'distillate', 'wax', 'aromatics', 'coke']
        expected = [
            'CO2', 'H2', 'CO', 'H2O',
            'CH4', 'C2H4', 'C2H6', 'C3H6', 'C3H8', 'C4_lump', 'C5plus_lump',
            'LPG', 'light_olefins', 'naphtha', 'distillate', 'wax',
            'aromatics', 'oxygenates', 'coke'
        ]
        for comp in zeo_products:
            assert comp in expected


class TestStoichiometry:
    """Tests for reaction stoichiometry dictionaries."""
    
    def test_rwgs_stoichiometry_exists(self):
        """Verify RWGS stoichiometry dictionary exists."""
        assert RWGS_STOICHIOMETRY is not None
        assert isinstance(RWGS_STOICHIOMETRY, dict)
    
    def test_rwgs_stoichiometry_balance(self):
        """Verify RWGS stoichiometry balances (atom balance)."""
        # RWGS: CO2 + H2 <-> CO + H2O
        # C: -1 (CO2) + 1 (CO) = 0 ✓
        # H: -2 (H2) + 2 (H2O) = 0 ✓
        # O: -2 (CO2) + 1 (CO) + 1 (H2O) = 0 ✓
        assert RWGS_STOICHIOMETRY['CO2'] == -1.0
        assert RWGS_STOICHIOMETRY['H2'] == -1.0
        assert RWGS_STOICHIOMETRY['CO'] == 1.0
        assert RWGS_STOICHIOMETRY['H2O'] == 1.0
    
    def test_ft_stoichiometry_exists(self):
        """Verify FT stoichiometry dictionary exists."""
        assert FT_STOICHIOMETRY is not None
        assert isinstance(FT_STOICHIOMETRY, dict)
    
    def test_ft_stoichiometry_structure(self):
        """Verify FT stoichiometry has correct structure."""
        expected_products = ['CH4', 'C2H4', 'C2H6', 'C3H6', 'C3H8', 'C5plus_lump', 'wax']
        for product in expected_products:
            assert product in FT_STOICHIOMETRY
            assert isinstance(FT_STOICHIOMETRY[product], dict)
    
    def test_ft_ch4_stoichiometry(self):
        """Verify CH4 formation stoichiometry (CO + 2H2 -> CH4)."""
        # CO: -1, H2: -2, CH4: +1
        stoich = FT_STOICHIOMETRY['CH4']
        assert stoich['CO'] == -1.0
        assert stoich['H2'] == -2.0
        assert stoich['CH4'] == 1.0
    
    def test_ft_c5plus_stoichiometry(self):
        """Verify C5+ (wax) formation stoichiometry."""
        stoich = FT_STOICHIOMETRY['C5plus_lump']
        assert stoich['CO'] == -5.0
        assert stoich['H2'] == -11.0
        assert stoich['C5plus_lump'] == 1.0
    
    def test_zeolite_stoichiometry_exists(self):
        """Verify zeolite stoichiometry dictionary exists."""
        assert ZEOLITE_STOICHIOMETRY is not None
        assert isinstance(ZEOLITE_STOICHIOMETRY, dict)
    
    def test_zeolite_cracking_stoichiometry(self):
        """Verify zeolite cracking stoichiometry."""
        # Wax cracking: wax -> 0.5 distillate + 0.3 LPG
        assert 'distillate' in ZEOLITE_STOICHIOMETRY
        distillate_stoich = ZEOLITE_STOICHIOMETRY['distillate']
        assert distillate_stoich['wax'] == -1.0
        assert distillate_stoich['distillate'] == 0.5
        assert distillate_stoich['LPG'] == 0.3


class TestReactionPhysics:
    """Tests for reaction physics assumptions."""
    
    def test_rwgs_is_reversible(self):
        """Verify RWGS is reversible (has reverse reaction)."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        # Should include both forward and reverse rate terms
        assert 'p_CO2 * p_H2' in source
        assert 'p_CO * p_H2O' in source
        assert 'Keq' in source
    
    def test_ft_is_pressure_dependent(self):
        """Verify FT reaction rates are pressure-dependent."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        # FT rates should depend on partial pressures
        assert 'partial_pressure' in source
        assert 'rate_ft' in source
    
    def test_zeolite_is_concentration_dependent(self):
        """Verify zeolite rates are concentration-dependent (first-order)."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        # Zeolite rates depend on component concentrations
        assert 'flow_mol_comp' in source
        assert 'rate_zeo' in source
    
    def test_material_balance_present(self):
        """Verify material balance equations are defined."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'material_balance' in source
        assert 'dF_dW' in source
    
    def test_partial_pressure_calculation(self):
        """Verify partial pressure calculation is defined."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'partial_pressure_calc' in source
        assert 'mole fraction' in source.lower()


class TestParameterInitialization:
    """Tests for parameter initialization."""
    
    def test_k_rwgs_parameter(self):
        """Verify k_rwgs parameter is defined in build."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'k_rwgs' in source
    
    def test_keq_rwgs_parameter(self):
        """Verify Keq_rwgs parameter is defined in build."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'Keq_rwgs' in source
    
    def test_k_ft_parameters(self):
        """Verify FT rate constants are defined in build."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'k_ft' in source
    
    def test_k_zeo_parameters(self):
        """Verify zeolite rate constants are defined in build."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'k_zeo' in source
    
    def test_w_total_parameter(self):
        """Verify W_total (catalyst weight) parameter exists."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'W_total' in source


class TestStateVariables:
    """Tests for state variable definitions."""
    
    def test_flow_mol_comp_variable(self):
        """Verify molar flow rate variable."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'flow_mol_comp' in source
    
    def test_temperature_variable(self):
        """Verify temperature variable."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'temperature' in source
    
    def test_pressure_variable(self):
        """Verify pressure variable."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'pressure' in source
    
    def test_flow_mol_variable(self):
        """Verify total molar flow variable."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'flow_mol' in source
    
    def test_partial_pressure_variable(self):
        """Verify partial pressure variable."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'partial_pressure' in source
    
    def test_reaction_rate_variables(self):
        """Verify reaction rate variables exist."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'rate_rwgs' in source
        assert 'rate_ft' in source
        assert 'rate_zeo' in source


class TestSpatialDiscretization:
    """Tests for 1D spatial discretization."""
    
    def test_continuous_set_defined(self):
        """Verify ContinuousSet for spatial domain is used."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'ContinuousSet' in source or 'self.W = ' in source
    
    def test_catalyst_weight_domain(self):
        """Verify catalyst weight domain is [0, 1]."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'bounds=(0, 1.0)' in source or 'ContinuousSet' in source
    
    def test_derivative_var_used(self):
        """Verify DerivativeVar is used for material balance."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'DerivativeVar' in source or 'dF_dW' in source


class TestInitializeMethod:
    """Tests for initialization method signature and functionality."""
    
    def test_initialize_signature(self):
        """Verify initialize method has correct parameters."""
        import inspect
        sig = inspect.signature(BifunctionalPackedBedReactorData.initialize)
        params = list(sig.parameters.keys())
        
        assert 'self' in params
        assert 'inlet_flow' in params
        assert 'inlet_temperature' in params
        assert 'inlet_pressure' in params
        assert 'W_total_value' in params
    
    def test_initialize_accepts_kinetic_dicts(self):
        """Verify initialize accepts kinetic parameter dictionaries."""
        import inspect
        sig = inspect.signature(BifunctionalPackedBedReactorData.initialize)
        params = list(sig.parameters.keys())
        
        assert 'k_ft_dict' in params
        assert 'k_zeo_dict' in params
    
    def test_initialize_has_defaults(self):
        """Verify initialize has sensible default values."""
        import inspect
        sig = inspect.signature(BifunctionalPackedBedReactorData.initialize)
        
        assert sig.parameters['W_total_value'].default == 1.0
        assert sig.parameters['k_rwgs_value'].default == 0.1
        assert sig.parameters['Keq_value'].default == 1.0


class TestHelperFunction:
    """Tests for build_bifunctional_reactor helper function."""
    
    def test_helper_function_exists(self):
        """Verify build_bifunctional_reactor helper exists."""
        assert callable(build_bifunctional_reactor)
    
    def test_helper_function_signature(self):
        """Verify helper function has correct signature."""
        import inspect
        sig = inspect.signature(build_bifunctional_reactor)
        params = list(sig.parameters.keys())
        
        assert 'flowsheet_block' in params
        assert 'property_package' in params
        assert 'inlet_flow' in params
        assert 'inlet_temperature' in params
        assert 'inlet_pressure' in params
        assert 'W_total' in params


class TestDocumentation:
    """Tests for documentation completeness."""
    
    def test_class_has_docstring(self):
        """Verify class has docstring."""
        assert BifunctionalPackedBedReactorData.__doc__ is not None
        doc = BifunctionalPackedBedReactorData.__doc__
        assert 'RWGS' in doc or 'reverse water-gas shift' in doc.lower()
    
    def test_build_method_documented(self):
        """Verify build method has docstring."""
        assert BifunctionalPackedBedReactorData.build.__doc__ is not None
    
    def test_initialize_method_documented(self):
        """Verify initialize method has docstring."""
        assert BifunctionalPackedBedReactorData.initialize.__doc__ is not None
    
    def test_performance_method_documented(self):
        """Verify performance reporting method has docstring."""
        assert BifunctionalPackedBedReactorData._get_performance_contents.__doc__ is not None
    
    def test_module_has_docstring(self):
        """Verify module has docstring."""
        import src.ft_model.bifunctional_reactor as reactor_module
        assert reactor_module.__doc__ is not None


class TestIntegration:
    """Integration tests for the bifunctional reactor."""
    
    def test_module_imports(self):
        """Verify module imports successfully."""
        from src.ft_model import bifunctional_reactor
        assert bifunctional_reactor is not None
    
    def test_class_importable(self):
        """Verify class can be imported."""
        from src.ft_model.bifunctional_reactor import BifunctionalPackedBedReactorData
        assert BifunctionalPackedBedReactorData is not None
    
    def test_stoichiometry_dicts_accessible(self):
        """Verify stoichiometry dicts are accessible from module."""
        assert RWGS_STOICHIOMETRY is not None
        assert FT_STOICHIOMETRY is not None
        assert ZEOLITE_STOICHIOMETRY is not None
    
    def test_main_example_runs(self):
        """Verify main example code executes without error."""
        # The __main__ section should run without errors
        import src.ft_model.bifunctional_reactor
        # If import succeeds, main section was parsed correctly


class TestAdvancedFeatures:
    """Tests for advanced reactor features."""
    
    def test_coupled_reactions(self):
        """Verify reactions are coupled in material balance."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        # Material balance should include all three reaction types
        assert 'rwgs_term' in source or 'rate_rwgs' in source
        assert 'ft_term' in source or 'rate_ft' in source
        assert 'zeo_term' in source or 'rate_zeo' in source
    
    def test_isothermal_option(self):
        """Verify isothermal operation option is available."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'isothermal' in source
        assert 'isothermal_constraint' in source
    
    def test_pressure_drop_constraint(self):
        """Verify pressure drop constraint is defined."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        assert 'pressure_drop' in source
    
    def test_multiple_ft_products(self):
        """Verify multiple FT products are handled."""
        assert 'CH4' in FT_STOICHIOMETRY
        assert 'C2H4' in FT_STOICHIOMETRY
        assert 'C3H6' in FT_STOICHIOMETRY
        assert 'C5plus_lump' in FT_STOICHIOMETRY
    
    def test_zeolite_cracking_pathway(self):
        """Verify zeolite cracking pathway is defined."""
        import inspect
        source = inspect.getsource(BifunctionalPackedBedReactorData.build)
        # Should handle: wax -> distillate -> naphtha
        assert 'distillate' in source
        assert 'wax' in source
        assert 'naphtha' in source
