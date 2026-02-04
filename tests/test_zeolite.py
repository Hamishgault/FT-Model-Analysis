"""
Unit Tests for Zeolite Models

Tests for zeolite reactor kinetics, unit model construction, and integration.
"""

import pytest
from src.zeolite_model import kinetics, reactor
from src.utils import lumping


class TestZeoliteKinetics:
    """Tests for zeolite kinetics module."""
    
    def test_import_kinetics(self):
        """Test that kinetics module imports successfully."""
        assert kinetics is not None
    
    def test_cracking_rate_exists(self):
        """Test that cracking rate function exists."""
        assert hasattr(kinetics, 'cracking_rate')
    
    def test_isomerization_rate_exists(self):
        """Test that isomerization rate function exists."""
        assert hasattr(kinetics, 'isomerization_rate')
    
    def test_oligomerization_rate_exists(self):
        """Test that oligomerization rate function exists."""
        assert hasattr(kinetics, 'oligomerization_rate')
    
    def test_deactivation_rate_exists(self):
        """Test that deactivation rate function exists."""
        assert hasattr(kinetics, 'deactivation_rate')
    
    def test_selectivity_function_exists(self):
        """Test that selectivity function exists."""
        assert hasattr(kinetics, 'selectivity_to_gasoline')


class TestZeoliteReactor:
    """Tests for zeolite reactor unit models."""
    
    def test_import_reactor(self):
        """Test that reactor module imports successfully."""
        assert reactor is not None
    
    def test_zeolite_reactor_class_exists(self):
        """Test that ZeoliteReactor class exists."""
        assert hasattr(reactor, 'ZeoliteReactor')
    
    def test_zeolite_reactor_data_class_exists(self):
        """Test that ZeoliteReactorData class exists."""
        assert hasattr(reactor, 'ZeoliteReactorData')
    
    def test_build_zeolite_reactor_function_exists(self):
        """Test that build_zeolite_reactor factory function exists."""
        assert hasattr(reactor, 'build_zeolite_reactor')
    
    def test_configure_kinetics_function_exists(self):
        """Test that configure_reactor_kinetics function exists."""
        assert hasattr(reactor, 'configure_reactor_kinetics')
    
    def test_heat_removal_function_exists(self):
        """Test that heat removal function exists."""
        assert hasattr(reactor, 'add_heat_removal')


class TestLumping:
    """Tests for lumping utilities."""
    
    def test_import_lumping(self):
        """Test that lumping module imports successfully."""
        assert lumping is not None
    
    def test_asf_to_lumped_function_exists(self):
        """Test that ASF to lumped function exists."""
        assert hasattr(lumping, 'asf_distribution_to_lumped')
    
    def test_ft_to_zeolite_feed_function_exists(self):
        """Test that FT to zeolite conversion function exists."""
        assert hasattr(lumping, 'ft_product_to_zeolite_feed')
    
    def test_component_groups_function_exists(self):
        """Test that component groups function exists."""
        assert hasattr(lumping, 'create_lumped_component_groups')
    
    def test_carbon_lumping_function_exists(self):
        """Test that carbon number lumping function exists."""
        assert hasattr(lumping, 'apply_carbon_number_lumping')


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
