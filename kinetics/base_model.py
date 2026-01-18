"""Abstract base class for FT kinetic models.

Defines a minimal interface that kinetic model implementations should follow.
"""
from abc import ABC, abstractmethod
from typing import Sequence
import numpy as np


class KineticModel(ABC):
    """Base class for kinetic models.

    Implementations must provide `rate(T, P, Fi)` which returns an array-like
    of reaction rates r_j [mol m^-3 s^-1] for each reaction j.
    """

    @abstractmethod
    def rate(self, T: float, P: float, Fi: np.ndarray) -> np.ndarray:
        """Compute reaction rates.

        Parameters
        ----------
        T : float
            Temperature in K
        P : float
            Pressure in Pa
        Fi : np.ndarray
            Current molar flows of species (shape: n_species, units: mol/s)

        Returns
        -------
        np.ndarray
            Reaction rates for each reaction (shape: n_rxns, units: mol m^-3 s^-1)
        """
        raise NotImplementedError
