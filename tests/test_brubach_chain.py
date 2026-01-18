import numpy as np

from kinetics.brubach_2022 import BrubachModel, BrubachParams
from utils.species import SPECIES_IDX


def test_brubach_chain_small_nmax():
    params = BrubachParams(n_max=6)
    model = BrubachModel(params)

    Fi = np.zeros(len(SPECIES_IDX))
    Fi[SPECIES_IDX["CO2"]] = 1.0
    Fi[SPECIES_IDX["H2"]] = 3.0

    # should solve surface and return rates
    r = model.rate(523.15, 1e5, Fi)
    r = np.asarray(r)
    assert r.shape[0] == model.n_rxns
    assert np.isfinite(r).all()

    # run PFR for a short reactor
    from reactor.pfr import PFR
    pfr = PFR(A=1e-3, L=0.1, model=model, nu=None)
    z = np.linspace(0.0, 0.1, 11)
    sol = pfr.run(523.15, 1e5, Fi, z_eval=z)
    # flows finite and non-negative
    assert np.all(np.isfinite(sol.y))
    assert np.all(sol.y >= -1e-12)
