"""Parameter sets and reactor defaults for example simulations."""
from typing import Dict
import numpy as np

from .species import SPECIES_IDX, NU

# Reactor defaults
REACTOR = {
    "A": 1e-3,  # m^2
    "L": 5.0,  # m
}

# Inlet molar flows (mol/s)
# CO:H2 typical syngas ratio ~1:3 for FT
F_INLET = np.zeros(len(SPECIES_IDX), dtype=float)
F_INLET[SPECIES_IDX["CO"]] = 1.0
F_INLET[SPECIES_IDX["H2"]] = 3.0

# Temperature and pressure
T0 = 523.15  # K (~250 C)
P0 = 1e5  # Pa

# Example model parameter sets (placeholders)
MODEL_PARAMS: Dict[str, dict] = {
    "power_law": {"k": 1e-3, "n_CO": 1.0, "n_H2": 0.5},
    "iglesia": {"k0": 1e-3, "k1": 5e-4},
    "steynberg": {"k": 2e-3, "K_CO": 1.0},
    "vanderlaan": {"k": 1.2e-3, "alpha_H2": 0.2},
}
