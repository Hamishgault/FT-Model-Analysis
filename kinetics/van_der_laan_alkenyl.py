"""Placeholder: Van der Laan alkenyl mechanism (simplified).

One effective reaction with a weak H2 dependence (placeholder).
"""



from typing import Optional
import numpy as np

from .base_model import KineticModel
from utils.species import SPECIES_IDX


class VanDerLaanAlkenylModel(KineticModel):
    def __init__(self, params: Optional[dict] = None):
        self.params = params or {"k": 1.2e-3, "alpha_H2": 0.2}
        self.n_rxns = 1

    def rate(self, T: float, P: float, Fi: np.ndarray) -> np.ndarray:
        iCO = SPECIES_IDX["CO"]
        iH2 = SPECIES_IDX["H2"]
        k = self.params["k"]
        alpha = self.params["alpha_H2"]

        CO = max(Fi[iCO], 0.0)
        H2 = max(Fi[iH2], 0.0)
        r = k * CO * (1.0 + alpha * H2)
        return np.array([r])