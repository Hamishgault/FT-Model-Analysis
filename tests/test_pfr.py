import numpy as np
from kinetics.base_model import KineticModel
from reactor.pfr import PFR
from utils.species import NU, SPECIES_IDX


class ZeroModel(KineticModel):
    def __init__(self):
        self.n_rxns = NU.shape[1]

    def rate(self, T, P, Fi):
        return np.zeros(self.n_rxns)


def test_pfr_with_zero_rates_returns_constant_flows():
    model = ZeroModel()
    pfr = PFR(A=1e-3, L=1.0, model=model, nu=NU)

    F0 = np.zeros(len(SPECIES_IDX))
    F0[SPECIES_IDX["CO"]] = 1.0
    F0[SPECIES_IDX["H2"]] = 3.0

    sol = pfr.run(T=300.0, P=1e5, F0=F0, z_eval=np.linspace(0.0, 1.0, 11))

    # with zero rates flows must remain unchanged
    assert sol.y.shape[0] == len(SPECIES_IDX)
    assert sol.y.shape[1] == 11
    assert np.allclose(sol.y[:, -1], F0)
