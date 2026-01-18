"""Placeholder: Steynberg carbide mechanism (high-level placeholder).

Provides a single effective rate representing carbide-mediated CO conversion.
"""
from typing import Optional
import numpy as np

from .base_model import KineticModel
from utils.species import SPECIES_IDX


class SteynbergCarbideModel(KineticModel):
    """Simplified carbide mechanism placeholder.

    The model uses a Langmuir-type form (placeholder) with CO adsorption
    limitation.
    """

    def __init__(self, params: Optional[dict] = None):
        self.params = params or {"k": 2e-3, "K_CO": 1.0}
        self.n_rxns = 1

    def rate(self, T: float, P: float, Fi: np.ndarray) -> np.ndarray:
        iCO = SPECIES_IDX["CO"]
        k = self.params["k"]
        K_CO = self.params["K_CO"]

        CO = max(Fi[iCO], 0.0)
        r = k * CO / (1.0 + K_CO * CO)
        return np.array([r])