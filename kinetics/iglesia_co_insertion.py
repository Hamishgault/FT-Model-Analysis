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
        # k0: activation rate constant, k1: chain growth constant
        # fractions determine product distribution of chain growth
        self.params = params or {"k0": 1e-3, "k1": 5e-4, "f_CH4": 0.33, "f_C2_4": 0.33}
        # now return three product-forming rates to match global NU which has 3 columns
        self.n_rxns = 3

        # explicitly provide a nu that maps the three rates to species changes
        from utils.species import NU
        # NU already has 3 columns in utils.species, so reference it directly
        self.nu = NU

    def rate(self, T: float, P: float, Fi: np.ndarray) -> np.ndarray:
        """Return three product-forming rates: CH4, C2-4, C5+.

        Mechanistic placeholders:
        - r_activation = k0 * CO
        - r_chain = k1 * r_activation
        - distribute r_chain into product channels using fractions
        """
        iCO = SPECIES_IDX["CO"]
        k0 = self.params["k0"]
        k1 = self.params["k1"]
        f_CH4 = self.params.get("f_CH4", 0.33)
        f_C2_4 = self.params.get("f_C2_4", 0.33)

        CO = max(Fi[iCO], 0.0)
        r_activation = k0 * CO
        r_chain = k1 * r_activation

        r_ch4 = f_CH4 * r_chain
        r_c2_4 = f_C2_4 * r_chain
        r_c5 = max(r_chain - r_ch4 - r_c2_4, 0.0)

        return np.array([r_ch4, r_c2_4, r_c5])
