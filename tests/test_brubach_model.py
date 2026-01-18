import numpy as np

from kinetics.brubach_2022 import BrubachModel
from reactor.pfr import PFR
from utils.species import SPECIES_IDX, NU


def test_brubach_rate_shape_and_non_negative():
    model = BrubachModel()
    Fi = np.zeros(len(SPECIES_IDX))
    Fi[SPECIES_IDX["CO"]] = 1.0
    Fi[SPECIES_IDX["H2"]] = 3.0
    Fi[SPECIES_IDX["CO2"]] = 0.1

    r = model.rate(523.15, 1e5, Fi)
    r = np.asarray(r)

    assert r.ndim == 1
    # product-forming reaction rates (indices 5,6,7) should be non-negative
    assert r[5] >= -1e-12 and r[6] >= -1e-12 and r[7] >= -1e-12
    # other rates may be negative (net desorption) but must be finite
    assert np.isfinite(r).all()


def test_brubach_runs_in_pfr():
    model = BrubachModel()
    pfr = PFR(A=1e-3, L=1.0, model=model, nu=None)

    F0 = np.zeros(len(SPECIES_IDX))
    F0[SPECIES_IDX["CO"]] = 1.0
    F0[SPECIES_IDX["H2"]] = 3.0
    F0[SPECIES_IDX["CO2"]] = 0.1

    z_eval = np.linspace(0.0, 1.0, 11)
    sol = pfr.run(T=523.15, P=1e5, F0=F0, z_eval=z_eval)

    assert sol.y.shape[0] == len(SPECIES_IDX)
    assert sol.y.shape[1] == len(z_eval)
    # basic sanity checks: non-negative flows and finite values
    assert np.all(sol.y >= -1e-12)
    assert np.isfinite(sol.y).all()
