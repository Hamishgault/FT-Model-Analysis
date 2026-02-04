"""
Component Data Model

Defines the Component class for representing thermodynamic properties
and molecular characteristics of chemical species.
"""

from typing import Dict, List, Optional


class Component:
    """
    Represents a chemical component with thermodynamic properties.
    
    Attributes
    ----------
    name : str
        Component identifier/name
    molecular_weight : float
        Molecular weight [kg/mol]
    critical_temperature : float
        Critical temperature [K]
    critical_pressure : float
        Critical pressure [Pa]
    acentric_factor : float
        Acentric factor [-]
    heat_capacity : dict
        Heat capacity coefficients {a, b, c, d} for polynomial expansion
        Cp = a + b*T + c*T^2 + d*T^3
    phase : str
        Standard phase at room conditions ('vap', 'liq', 'solid')
    """
    
    REQUIRED_FIELDS = [
        'molecular_weight',
        'critical_temperature',
        'critical_pressure',
        'acentric_factor',
        'phase'
    ]
    
    VALID_PHASES = {'vap', 'liq', 'solid'}
    
    def __init__(
        self,
        name: str,
        molecular_weight: float,
        critical_temperature: float,
        critical_pressure: float,
        acentric_factor: float,
        phase: str,
        heat_capacity: Optional[Dict[str, float]] = None,
        **kwargs
    ):
        """
        Initialize a Component.
        
        Parameters
        ----------
        name : str
            Component name
        molecular_weight : float
            Molecular weight [kg/mol]
        critical_temperature : float
            Critical temperature [K]
        critical_pressure : float
            Critical pressure [Pa]
        acentric_factor : float
            Acentric factor [-]
        phase : str
            Phase ('vap', 'liq', 'solid')
        heat_capacity : dict, optional
            Heat capacity coefficients
        **kwargs
            Additional properties
        """
        self.name = name
        self.molecular_weight = molecular_weight
        self.critical_temperature = critical_temperature
        self.critical_pressure = critical_pressure
        self.acentric_factor = acentric_factor
        self.phase = phase
        self.heat_capacity = heat_capacity or {}
        self.additional_properties = kwargs
        
        # Validate on initialization
        self.validate()
    
    def validate(self) -> None:
        """
        Validate that all required fields are present and valid.
        
        Raises
        ------
        ValueError
            If required fields are missing or invalid
        """
        # Check required fields
        if self.phase not in self.VALID_PHASES:
            raise ValueError(
                f"Invalid phase '{self.phase}'. Must be one of {self.VALID_PHASES}"
            )
        
        # Check numeric fields are positive
        if self.molecular_weight <= 0:
            raise ValueError(f"Molecular weight must be positive, got {self.molecular_weight}")
        
        if self.critical_temperature <= 0:
            raise ValueError(f"Critical temperature must be positive, got {self.critical_temperature}")
        
        if self.critical_pressure <= 0:
            raise ValueError(f"Critical pressure must be positive, got {self.critical_pressure}")
    
    def to_dict(self) -> Dict:
        """
        Convert component to dictionary representation.
        
        Returns
        -------
        dict
            Component data as dictionary
        """
        return {
            'name': self.name,
            'molecular_weight': self.molecular_weight,
            'critical_temperature': self.critical_temperature,
            'critical_pressure': self.critical_pressure,
            'acentric_factor': self.acentric_factor,
            'phase': self.phase,
            'heat_capacity': self.heat_capacity,
            **self.additional_properties
        }
    
    def to_idaes_format(self) -> Dict:
        """
        Convert component to IDAES Generic Property Package format.
        
        Returns
        -------
        dict
            Component data formatted for IDAES
        """
        return {
            'molecular_weight': (self.molecular_weight, 'kg/mol'),
            'critical_temperature': (self.critical_temperature, 'K'),
            'critical_pressure': (self.critical_pressure, 'Pa'),
            'acentric_factor': (self.acentric_factor, None),
        }
    
    def __repr__(self) -> str:
        return f"Component(name='{self.name}', MW={self.molecular_weight:.3f}, phase={self.phase})"
    
    def __str__(self) -> str:
        return self.name
