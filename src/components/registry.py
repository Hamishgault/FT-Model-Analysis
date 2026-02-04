"""
Component Registry

Manages loading and access to component definitions from YAML files.
Provides a centralized registry for all chemical species used in the process model.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional

import yaml

from .component import Component


class ComponentRegistry:
    """
    Registry for managing chemical components loaded from YAML files.
    
    Loads YAML files from the data/ directory and provides methods to:
    - Retrieve components by name
    - List all components or filter by phase
    - Export to IDAES-compatible format
    """
    
    def __init__(self, data_dir: Optional[str] = None):
        """
        Initialize the component registry.
        
        Parameters
        ----------
        data_dir : str, optional
            Path to directory containing YAML component files.
            If None, uses src/components/data/
        """
        if data_dir is None:
            # Get the data directory relative to this file
            current_dir = Path(__file__).parent
            data_dir_path: Path = current_dir / 'data'
        else:
            data_dir_path = Path(data_dir)
        
        self.data_dir: Path = data_dir_path
        self._components: Dict[str, Component] = {}
        self._load_all_components()
    
    def _load_all_components(self) -> None:
        """Load all YAML files from the data directory."""
        if not self.data_dir.exists():
            raise FileNotFoundError(f"Component data directory not found: {self.data_dir}")
        
        yaml_files = sorted(self.data_dir.glob('*.yaml'))
        
        for yaml_file in yaml_files:
            try:
                self._load_yaml_file(yaml_file)
            except Exception as e:
                raise RuntimeError(f"Error loading {yaml_file.name}: {str(e)}")
    
    def _load_yaml_file(self, yaml_file: Path) -> None:
        """
        Load components from a YAML file.
        
        Parameters
        ----------
        yaml_file : Path
            Path to YAML file
        """
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
        
        if data is None:
            return
        
        for name, props in data.items():
            # Ensure name is set
            props_copy = props.copy() if isinstance(props, dict) else {}
            if 'name' not in props_copy:
                props_copy['name'] = name
            
            try:
                component = Component(**props_copy)
                self._components[name] = component
            except Exception as e:
                raise ValueError(f"Error creating component '{name}': {str(e)}")
    
    def get(self, name: str) -> Optional[Component]:
        """
        Get a component by name.
        
        Parameters
        ----------
        name : str
            Component name
        
        Returns
        -------
        Component or None
            Component object if found, None otherwise
        """
        return self._components.get(name)
    
    def get_or_raise(self, name: str) -> Component:
        """
        Get a component by name, raising exception if not found.
        
        Parameters
        ----------
        name : str
            Component name
        
        Returns
        -------
        Component
            Component object
        
        Raises
        ------
        KeyError
            If component not found
        """
        if name not in self._components:
            available = ', '.join(sorted(self._components.keys()))
            raise KeyError(
                f"Component '{name}' not found. Available: {available}"
            )
        return self._components[name]
    
    def list_all(self) -> List[str]:
        """
        List all registered component names.
        
        Returns
        -------
        list
            Sorted list of component names
        """
        return sorted(self._components.keys())
    
    def list_by_phase(self, phase: str) -> List[str]:
        """
        List components by phase.
        
        Parameters
        ----------
        phase : str
            Phase filter ('vap', 'liq', 'solid')
        
        Returns
        -------
        list
            Component names with specified phase, sorted
        """
        return sorted([
            name for name, comp in self._components.items()
            if comp.phase == phase
        ])
    
    def to_dict(self) -> Dict[str, Dict]:
        """
        Convert all components to dictionary format.
        
        Returns
        -------
        dict
            All components as {name: component_dict}
        """
        return {name: comp.to_dict() for name, comp in self._components.items()}
    
    def to_idaes_property_package(self) -> Dict[str, Dict]:
        """
        Convert all components to IDAES Generic Property Package format.
        
        This creates a dictionary structure compatible with IDAES property
        package definitions.
        
        Returns
        -------
        dict
            Components formatted for IDAES {name: idaes_format_dict}
        """
        return {
            name: comp.to_idaes_format()
            for name, comp in self._components.items()
        }
    
    def get_component_count(self) -> int:
        """
        Get total number of registered components.
        
        Returns
        -------
        int
            Number of components
        """
        return len(self._components)
    
    def __repr__(self) -> str:
        count = self.get_component_count()
        return f"ComponentRegistry({count} components from {self.data_dir.name})"
