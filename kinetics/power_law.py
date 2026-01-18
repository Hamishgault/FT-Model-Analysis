"""Simple empirical power-law kinetic model (placeholder).

This file contains a minimal `PowerLawModel` implementation that follows the
`KineticModel` interface. The formulas are placeholders to allow the PFR solver
and comparison machinery to exercise the interface.
"""
from typing import Optional
import numpy as np

from .base_model import KineticModel
from utils.species import SPECIES_IDX


class PowerLawModel(KineticModel):
    """Empirical power-law model.

    Reaction scheme: single effective FT reaction consuming CO and H2 and
    producing a lumped hydrocarbon product (placeholder).
    """

    def __init__(self, params: Optional[dict] = None):
        # params: k [mol m^-3 s^-1 / (mol/s)^n], n_CO, n_H2
        self.params = params or {"k": 1e-3, "n_CO": 1.0, "n_H2": 0.5}
        self.n_rxns = 1

    def rate(self, T: float, P: float, Fi: np.ndarray) -> np.ndarray:
        """Compute single effective rate (placeholder).

        Notes
        -----
        - For now treat `Fi` (mol/s) as a proxy for concentration scaling. In a
          later iteration we'd divide by volumetric flow to get concentration.
        """
        k = self.params["k"]
        n_CO = self.params["n_CO"]
        n_H2 = self.params["n_H2"]

        iCO = SPECIES_IDX["CO"]
        iH2 = SPECIES_IDX["H2"]

        # simple empirical form
        r = k * (max(Fi[iCO], 0.0) ** n_CO) * (max(Fi[iH2], 0.0) ** n_H2)
        return np.array([r])
