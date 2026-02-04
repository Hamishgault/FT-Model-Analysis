"""
Unit Tests for Fischer-Tropsch Models

Tests for FT reactor kinetics, unit model construction, and integration.
"""

import pytest
from src.ft_model import kinetics, reactor
from src.utils import parameters


class TestFTKinetics:
    """Tests for FT kinetics module."""
    
    def test_import_kinetics(self):
        """Test that kinetics module imports successfully."""
        assert kinetics is not None
    
    def test_ft_rate_expression_exists(self):
        """Test that FT rate expression function exists."""
        assert hasattr(kinetics, 'ft_rate_expression')
    
    def test_co_consumption_rate_exists(self):
        """Test that CO consumption rate function exists."""
        assert hasattr(kinetics, 'co_consumption_rate')
    
    def test_h2_consumption_rate_exists(self):
        """Test that H2 consumption rate function exists."""
        assert hasattr(kinetics, 'h2_consumption_rate')
    
    def test_product_distribution_exists(self):
        """Test that product distribution function exists."""
        assert hasattr(kinetics, 'hydrocarbon_product_distribution')
    
    def test_selectivity_function_exists(self):
        """Test that selectivity function exists."""
        assert hasattr(kinetics, 'selectivity_to_liquid_hydrocarbons')


class TestFTReactor:
    """Tests for FT reactor unit models."""
    
    def test_import_reactor(self):
        """Test that reactor module imports successfully."""
        assert reactor is not None
    
    def test_ft_reactor_class_exists(self):
        """Test that FTReactor class exists."""
        assert hasattr(reactor, 'FTReactor')
    
    def test_ft_reactor_data_class_exists(self):
        """Test that FTReactorData class exists."""
        assert hasattr(reactor, 'FTReactorData')
    
    def test_build_ft_reactor_function_exists(self):
        """Test that build_ft_reactor factory function exists."""
        assert hasattr(reactor, 'build_ft_reactor')
    
    def test_configure_kinetics_function_exists(self):
        """Test that configure_reactor_kinetics function exists."""
        assert hasattr(reactor, 'configure_reactor_kinetics')
    
    def test_energy_balance_function_exists(self):
        """Test that energy balance function exists."""
        assert hasattr(reactor, 'add_energy_balance')


class TestIntegration:
    """Integration tests for FT models."""
    
    def test_parameters_module_imports(self):
        """Test that parameters module imports successfully."""
        assert parameters is not None
    
    def test_load_parameters_function_exists(self):
        """Test that load_parameters function exists."""
        assert hasattr(parameters, 'load_parameters')


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
