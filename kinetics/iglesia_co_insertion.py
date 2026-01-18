"""Placeholder: Iglesia's CO insertion style model (simplified).

This implementation is intentionally minimal and returns a vector of rates for
an illustrative two-step mechanism: (1) CO activation, (2) chain growth
(lumped). Both rates are placeholders.
"""
from typing import Optional
import numpy as np

from .base_model import KineticModel
from utils.species import SPECIES_IDX


class IglesiaCOInsertionModel(KineticModel):
    """Simplified two-step CO-insertion model (placeholder).

    - r0: CO activation (consumes CO)
    - r1: chain growth (consumes surface intermediates; implemented as
          proportional to CO activation for simplicity)
    """

    def __init__(self, params: Optional[dict] = None):
        self.params = params or {"k0": 1e-3, "k1": 5e-4}
        self.n_rxns = 2

    def rate(self, T: float, P: float, Fi: np.ndarray) -> np.ndarray:
        iCO = SPECIES_IDX["CO"]
        k0 = self.params["k0"]
        k1 = self.params["k1"]

        r0 = k0 * max(Fi[iCO], 0.0)
        r1 = k1 * r0  # placeholder coupling
        return np.array([r0, r1])
