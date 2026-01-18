import numpy as np
import pytest

from kinetics import (
    PowerLawModel,
    IglesiaCOInsertionModel,
    SteynbergCarbideModel,
    VanDerLaanAlkenylModel,
)
from utils.species import SPECIES_IDX


@pytest.mark.parametrize("model_cls", [
    PowerLawModel,
    IglesiaCOInsertionModel,
    SteynbergCarbideModel,
    VanDerLaanAlkenylModel,
])
def test_model_rate_shape_and_non_negative(model_cls):
    params = None
    model = model_cls(params)

    Fi = np.zeros(len(SPECIES_IDX))
    Fi[SPECIES_IDX["CO"]] = 1.0
    Fi[SPECIES_IDX["H2"]] = 3.0

    r = model.rate(523.15, 1e5, Fi)
    r = np.asarray(r)

    assert r.ndim == 1
    assert (r >= 0).all()
