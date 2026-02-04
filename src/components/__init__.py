"""
Components Module

Provides access to the component database system including Component class
and ComponentRegistry for managing thermodynamic properties.
"""

from .component import Component
from .registry import ComponentRegistry

__all__ = ['Component', 'ComponentRegistry']
