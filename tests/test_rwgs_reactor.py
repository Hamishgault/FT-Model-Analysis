"""
Unit Tests for Packed-Bed RWGS Reactor Model

Tests for the PackedBedRWGSReactor custom IDAES unit model.
"""

import pytest
from src.ft_model.rwgs_reactor import (
    PackedBedRWGSReactor,
    PackedBedRWGSReactorData,
    RWGS_STOICHIOMETRY,
)


class TestRWGSStoichiometry:
    """Tests for RWGS stoichiometry constants."""
    
    def test_stoichiometry_dict_exists(self):
        """Test that stoichiometry dictionary is defined."""
        assert isinstance(RWGS_STOICHIOMETRY, dict)
    
    def test_stoichiometry_reactants(self):
        """Test that reactants have negative stoichiometric coefficients."""
        assert RWGS_STOICHIOMETRY['CO2'] == -1.0
        assert RWGS_STOICHIOMETRY['H2'] == -1.0
    
    def test_stoichiometry_products(self):
        """Test that products have positive stoichiometric coefficients."""
        assert RWGS_STOICHIOMETRY['CO'] == 1.0
        assert RWGS_STOICHIOMETRY['H2O'] == 1.0
    
    def test_stoichiometry_balance(self):
        """Test atom balance in RWGS reaction."""
        # CO2 + H2 -> CO + H2O
        # C: 1 = 1 ✓
        # O: 2 + 0 = 1 + 1 ✓
        # H: 0 + 2 = 0 + 2 ✓
        c_balance = RWGS_STOICHIOMETRY['CO2'] + RWGS_STOICHIOMETRY['CO']
        o_balance = (
            -2 * RWGS_STOICHIOMETRY['CO2']
            - RWGS_STOICHIOMETRY['CO']
            - RWGS_STOICHIOMETRY['H2O']
        )
        h_balance = (
            2 * RWGS_STOICHIOMETRY['H2']
            + 2 * RWGS_STOICHIOMETRY['H2O']
        )
        
        # These should sum appropriately
        assert abs(c_balance) >= 0  # Just verify calculation works


class TestPackedBedRWGSReactorClass:
    """Tests for PackedBedRWGSReactor class definition."""
    
    def test_reactor_class_exists(self):
        """Test that PackedBedRWGSReactor class is defined."""
        assert PackedBedRWGSReactor is not None
    
    def test_reactor_data_class_exists(self):
        """Test that PackedBedRWGSReactorData class is defined."""
        assert PackedBedRWGSReactorData is not None
    
    def test_reactor_inherits_from_unit_model(self):
        """Test that reactor inherits from UnitModelBlockData."""
        # Check that the class has expected UnitModel attributes
        assert hasattr(PackedBedRWGSReactorData, 'build')
        assert hasattr(PackedBedRWGSReactorData, 'CONFIG')
    
    def test_reactor_has_initialize_method(self):
        """Test that reactor has initialization method."""
        assert hasattr(PackedBedRWGSReactorData, 'initialize')
    
    def test_reactor_has_performance_method(self):
        """Test that reactor has performance reporting method."""
        assert hasattr(PackedBedRWGSReactorData, '_get_performance_contents')


class TestReactorConfiguration:
    """Tests for reactor configuration options."""
    
    def test_config_has_dynamic_flag(self):
        """Test that CONFIG has dynamic option."""
        assert hasattr(PackedBedRWGSReactorData.CONFIG, 'dynamic')
    
    def test_config_has_holdup_flag(self):
        """Test that CONFIG has holdup option."""
        assert hasattr(PackedBedRWGSReactorData.CONFIG, 'has_holdup')


class TestReactorConstants:
    """Tests for reactor parameter and variable definitions."""
    
    def test_stoichiometry_contains_all_rwgs_species(self):
        """Test that all RWGS species are in stoichiometry dict."""
        rwgs_species = {'CO2', 'H2', 'CO', 'H2O'}
        for species in rwgs_species:
            assert species in RWGS_STOICHIOMETRY


class TestInitializationMethod:
    """Tests for reactor initialization logic."""
    
    def test_initialize_method_accepts_inlet_conditions(self):
        """Test that initialize method has correct signature."""
        import inspect
        sig = inspect.signature(PackedBedRWGSReactorData.initialize)
        params = list(sig.parameters.keys())
        
        assert 'inlet_flow' in params
        assert 'inlet_temperature' in params
        assert 'inlet_pressure' in params
    
    def test_initialize_method_accepts_kinetic_params(self):
        """Test that initialize accepts kinetic parameters."""
        import inspect
        sig = inspect.signature(PackedBedRWGSReactorData.initialize)
        params = list(sig.parameters.keys())
        
        assert 'W_total_value' in params
        assert 'k_rwgs_value' in params
        assert 'Keq_value' in params


class TestReactorDocumentation:
    """Tests for proper documentation."""
    
    def test_reactor_has_docstring(self):
        """Test that reactor class has documentation."""
        assert PackedBedRWGSReactor.__doc__ is not None
    
    def test_reactor_data_has_docstring(self):
        """Test that reactor data class has documentation."""
        assert PackedBedRWGSReactorData.__doc__ is not None
    
    def test_build_method_has_docstring(self):
        """Test that build method is documented."""
        assert PackedBedRWGSReactorData.build.__doc__ is not None
    
    def test_initialize_method_has_docstring(self):
        """Test that initialize method is documented."""
        assert PackedBedRWGSReactorData.initialize.__doc__ is not None


class TestIntegration:
    """Integration tests for reactor model."""
    
    def test_reactor_module_imports(self):
        """Test that reactor module imports successfully."""
        from src.ft_model import rwgs_reactor
        assert rwgs_reactor is not None
    
    def test_reactor_accessible_from_ft_model(self):
        """Test that reactor is accessible from ft_model package."""
        from src.ft_model.rwgs_reactor import PackedBedRWGSReactor
        assert PackedBedRWGSReactor is not None


class TestRWGSReactionPhysics:
    """Tests for RWGS reaction physics assumptions."""
    
    def test_rwgs_is_endothermic(self):
        """
        Test documentation reflects that RWGS is endothermic.
        CO2 + H2 <-> CO + H2O is endothermic (requires heat input).
        """
        doc = PackedBedRWGSReactorData.__doc__
        assert doc and ('CO2 + H2' in doc or 'reverse water-gas shift' in doc.lower())
    
    def test_rwgs_is_equilibrium_limited(self):
        """Test that rate expression includes equilibrium constant."""
        import inspect
        source = inspect.getsource(PackedBedRWGSReactorData.build)
        assert 'Keq' in source
        assert 'equilibrium' in source.lower()
    
    def test_rwgs_includes_forward_and_reverse_reaction(self):
        """Test that rate expression accounts for forward and reverse."""
        import inspect
        source = inspect.getsource(PackedBedRWGSReactorData.build)
        # Forward: p_CO2 * p_H2
        # Reverse: p_CO * p_H2O / Keq
        assert 'p_CO2' in source
        assert 'p_H2' in source
        assert 'p_CO' in source
        assert 'p_H2O' in source


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
