"""
Unit Tests for Component Database System

Tests for Component class, ComponentRegistry, and YAML loading.
"""

import pytest
from src.components import Component, ComponentRegistry


class TestComponent:
    """Tests for Component class."""
    
    def test_component_creation(self):
        """Test basic component creation."""
        comp = Component(
            name='CO2',
            molecular_weight=44.01e-3,
            critical_temperature=304.2,
            critical_pressure=7.38e6,
            acentric_factor=0.225,
            phase='vap'
        )
        assert comp.name == 'CO2'
        assert comp.molecular_weight == 44.01e-3
        assert comp.phase == 'vap'
    
    def test_component_with_heat_capacity(self):
        """Test component with heat capacity coefficients."""
        hc = {'a': 22.26, 'b': 5.981e-2, 'c': -3.501e-5, 'd': 7.469e-9}
        comp = Component(
            name='CO2',
            molecular_weight=44.01e-3,
            critical_temperature=304.2,
            critical_pressure=7.38e6,
            acentric_factor=0.225,
            phase='vap',
            heat_capacity=hc
        )
        assert comp.heat_capacity == hc
    
    def test_component_validation_invalid_phase(self):
        """Test that invalid phase raises error."""
        with pytest.raises(ValueError, match="Invalid phase"):
            Component(
                name='TestComp',
                molecular_weight=50.0e-3,
                critical_temperature=400.0,
                critical_pressure=5.0e6,
                acentric_factor=0.2,
                phase='invalid_phase'
            )
    
    def test_component_validation_negative_mw(self):
        """Test that negative molecular weight raises error."""
        with pytest.raises(ValueError, match="Molecular weight must be positive"):
            Component(
                name='TestComp',
                molecular_weight=-50.0e-3,
                critical_temperature=400.0,
                critical_pressure=5.0e6,
                acentric_factor=0.2,
                phase='vap'
            )
    
    def test_component_validation_negative_temp(self):
        """Test that negative critical temperature raises error."""
        with pytest.raises(ValueError, match="Critical temperature must be positive"):
            Component(
                name='TestComp',
                molecular_weight=50.0e-3,
                critical_temperature=-400.0,
                critical_pressure=5.0e6,
                acentric_factor=0.2,
                phase='vap'
            )
    
    def test_component_to_dict(self):
        """Test conversion to dictionary."""
        comp = Component(
            name='CO2',
            molecular_weight=44.01e-3,
            critical_temperature=304.2,
            critical_pressure=7.38e6,
            acentric_factor=0.225,
            phase='vap'
        )
        comp_dict = comp.to_dict()
        assert comp_dict['name'] == 'CO2'
        assert comp_dict['molecular_weight'] == 44.01e-3
        assert comp_dict['phase'] == 'vap'
    
    def test_component_to_idaes_format(self):
        """Test conversion to IDAES format."""
        comp = Component(
            name='CO2',
            molecular_weight=44.01e-3,
            critical_temperature=304.2,
            critical_pressure=7.38e6,
            acentric_factor=0.225,
            phase='vap'
        )
        idaes_format = comp.to_idaes_format()
        assert 'molecular_weight' in idaes_format
        assert idaes_format['molecular_weight'][0] == 44.01e-3
        assert idaes_format['molecular_weight'][1] == 'kg/mol'
    
    def test_component_repr(self):
        """Test string representation."""
        comp = Component(
            name='CO2',
            molecular_weight=44.01e-3,
            critical_temperature=304.2,
            critical_pressure=7.38e6,
            acentric_factor=0.225,
            phase='vap'
        )
        assert 'CO2' in repr(comp)
        assert 'Component' in repr(comp)


class TestComponentRegistry:
    """Tests for ComponentRegistry class."""
    
    @pytest.fixture
    def registry(self):
        """Create a registry instance for testing."""
        return ComponentRegistry()
    
    def test_registry_loads_components(self, registry):
        """Test that registry loads components from YAML files."""
        assert registry.get_component_count() > 0
    
    def test_registry_get_existing_component(self, registry):
        """Test retrieving an existing component."""
        co2 = registry.get('CO2')
        assert co2 is not None
        assert co2.name == 'CO2'
        assert co2.molecular_weight == 44.01e-3
    
    def test_registry_get_nonexistent_component(self, registry):
        """Test retrieving a non-existent component returns None."""
        comp = registry.get('NONEXISTENT')
        assert comp is None
    
    def test_registry_get_or_raise_existing(self, registry):
        """Test get_or_raise with existing component."""
        co2 = registry.get_or_raise('CO2')
        assert co2.name == 'CO2'
    
    def test_registry_get_or_raise_nonexistent(self, registry):
        """Test get_or_raise with non-existent component raises error."""
        with pytest.raises(KeyError):
            registry.get_or_raise('NONEXISTENT')
    
    def test_registry_list_all(self, registry):
        """Test listing all components."""
        all_comps = registry.list_all()
        assert isinstance(all_comps, list)
        assert len(all_comps) > 0
        assert 'CO2' in all_comps
    
    def test_registry_list_by_phase_vapor(self, registry):
        """Test listing components by phase (vapor)."""
        vapor_comps = registry.list_by_phase('vap')
        assert isinstance(vapor_comps, list)
        assert 'H2' in vapor_comps
        assert 'CO' in vapor_comps
    
    def test_registry_list_by_phase_liquid(self, registry):
        """Test listing components by phase (liquid)."""
        liquid_comps = registry.list_by_phase('liq')
        assert isinstance(liquid_comps, list)
        assert 'H2O' in liquid_comps
    
    def test_registry_to_dict(self, registry):
        """Test conversion to dictionary."""
        all_dict = registry.to_dict()
        assert isinstance(all_dict, dict)
        assert 'CO2' in all_dict
        assert all_dict['CO2']['name'] == 'CO2'
    
    def test_registry_to_idaes_property_package(self, registry):
        """Test conversion to IDAES property package format."""
        idaes_dict = registry.to_idaes_property_package()
        assert isinstance(idaes_dict, dict)
        assert 'CO2' in idaes_dict
        
        # Check IDAES format
        co2_idaes = idaes_dict['CO2']
        assert 'molecular_weight' in co2_idaes
        assert 'critical_temperature' in co2_idaes
        assert 'critical_pressure' in co2_idaes
        assert 'acentric_factor' in co2_idaes
    
    def test_registry_component_count(self, registry):
        """Test getting component count."""
        count = registry.get_component_count()
        assert isinstance(count, int)
        assert count > 0
    
    def test_registry_repr(self, registry):
        """Test registry string representation."""
        repr_str = repr(registry)
        assert 'ComponentRegistry' in repr_str
        assert 'components' in repr_str
    
    def test_registry_required_fields_co2(self, registry):
        """Test that required fields are present in CO2."""
        co2 = registry.get_or_raise('CO2')
        assert co2.molecular_weight is not None
        assert co2.critical_temperature is not None
        assert co2.critical_pressure is not None
        assert co2.acentric_factor is not None
        assert co2.phase is not None
    
    def test_registry_ft_components_exist(self, registry):
        """Test that FT synthesis components are loaded."""
        ft_components = ['CH4', 'C2H6', 'C2H4', 'C3H8', 'C3H6']
        for comp_name in ft_components:
            assert registry.get(comp_name) is not None, f"{comp_name} not found"
    
    def test_registry_zeolite_lumps_exist(self, registry):
        """Test that zeolite lumps are loaded."""
        zeolite_lumps = ['LPG', 'naphtha', 'distillate', 'wax', 'aromatics']
        for lump_name in zeolite_lumps:
            assert registry.get(lump_name) is not None, f"{lump_name} not found"


class TestComponentIntegration:
    """Integration tests for component system."""
    
    def test_full_workflow(self):
        """Test complete workflow from loading to IDAES export."""
        # Load registry
        registry = ComponentRegistry()
        
        # Get a component
        co2 = registry.get_or_raise('CO2')
        
        # Verify basic properties
        assert co2.name == 'CO2'
        assert co2.phase == 'vap'
        
        # Convert to IDAES format
        idaes_data = co2.to_idaes_format()
        assert isinstance(idaes_data, dict)
        
        # Export full property package
        property_pkg = registry.to_idaes_property_package()
        assert isinstance(property_pkg, dict)
        assert len(property_pkg) > 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
