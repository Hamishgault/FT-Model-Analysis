"""Species indices and stoichiometric matrix (placeholder).

Species:
    0: CO
    1: H2
    2: CH4
    3: C2_4 (lumped)
    4: C5plus (lumped heavy products)
    5: H2O

Reactions (placeholder):
    r0: formation of CH4: CO + 3 H2 -> CH4 + H2O
    r1: formation of C2-4: 2 CO + 5 H2 -> C2_4 + 2 H2O
    r2: formation of C5+: 5 CO + 12 H2 -> C5+ + 5 H2O

The stoichiometric matrix NU has shape (n_species, n_rxns) and entries
are nu_ij (change in species i per mole of reaction j).
"""
import numpy as np

SPECIES = ["CO", "H2", "CH4", "C2_4", "C5plus", "H2O", "CO2"]
SPECIES_IDX = {name: i for i, name in enumerate(SPECIES)}

# nu: rows = species, cols = reactions (placeholder global NU)
NU = np.array([
    # r0   r1   r2
    [-1.0, -2.0, -5.0],  # CO
    [-3.0, -5.0,-12.0],  # H2
    [1.0,   0.0,  0.0],  # CH4
    [0.0,   1.0,  0.0],  # C2_4
    [0.0,   0.0,  1.0],  # C5plus
    [1.0,   2.0,  5.0],  # H2O
    [0.0,   0.0,  0.0],  # CO2 (placeholder)
])


def compute_selectivity(F_in, F_out):
    """Compute product selectivities based on species molar flows.

    Selectivity defined as fraction of carbon converted to each product group
    from CO consumption. This is a placeholder and assumes each product
    contains the same number of carbon atoms per mole for normalization.
    """
    iCO = SPECIES_IDX["CO"]
    iCH4 = SPECIES_IDX["CH4"]
    iC2_4 = SPECIES_IDX["C2_4"]
    iC5 = SPECIES_IDX["C5plus"]

    CO_consumed = max(F_in[iCO] - F_out[iCO], 1e-12)

    prod_CH4 = max(F_out[iCH4] - F_in[iCH4], 0.0)
    prod_C2_4 = max(F_out[iC2_4] - F_in[iC2_4], 0.0)
    prod_C5 = max(F_out[iC5] - F_in[iC5], 0.0)

    total = prod_CH4 + prod_C2_4 + prod_C5
    if total <= 0:
        return {"CH4": 0.0, "C2_4": 0.0, "C5+": 0.0}

    return {
        "CH4": prod_CH4 / total,
        "C2_4": prod_C2_4 / total,
        "C5+": prod_C5 / total,
    }
