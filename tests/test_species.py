import numpy as np
import pytest

from utils.species import SPECIES_IDX, NU, compute_selectivity


def test_species_indices_and_nu_shape():
    assert len(SPECIES_IDX) == 7
    assert NU.shape == (7, 3)


def test_compute_selectivity_basic():
    F_in = np.zeros(len(SPECIES_IDX))
    F_in[SPECIES_IDX["CO"]] = 10.0

    F_out = F_in.copy()
    # produce 3 mol CH4 at the outlet
    F_out[SPECIES_IDX["CH4"]] += 3.0
    F_out[SPECIES_IDX["CO"]] -= 3.0

    sel = compute_selectivity(F_in, F_out)

    # all product carbon goes to CH4 in this simple case
    assert pytest.approx(1.0, rel=1e-6) == sel["CH4"]
    assert sel["C2_4"] == 0.0
    assert sel["C5+"] == 0.0
