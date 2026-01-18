import numpy as np
from reactor.pfr import PFR
from kinetics.base_model import KineticModel


class SingleRateModel(KineticModel):
    def __init__(self):
        # one reaction that consumes only CO
        self.nu = np.zeros((6, 1))
        from utils.species import SPECIES_IDX
        self.nu[SPECIES_IDX["CO"], 0] = -1.0
        self.nu[SPECIES_IDX["CH4"], 0] = 1.0

    def rate(self, T, P, Fi):
        return np.array([1.0])


def test_pfr_uses_model_nu_when_present():
    model = SingleRateModel()
    pfr = PFR(A=1e-3, L=1.0, model=model, nu=None)

    F0 = np.zeros(6)
    from utils.species import SPECIES_IDX
    F0[SPECIES_IDX["CO"]] = 2.0

    z, y = pfr.run(T=300.0, P=1e5, F0=F0, z_eval=np.linspace(0.0, 1.0, 3)).t, pfr.run(T=300.0, P=1e5, F0=F0, z_eval=np.linspace(0.0, 1.0, 3)).y

    # With constant rate 1.0 and nu mapping, CO should decrease along reactor
    assert y.shape[0] == 6
    assert y[SPECIES_IDX["CO"], -1] < F0[SPECIES_IDX["CO"]]
    assert y[SPECIES_IDX["CH4"], -1] > 0.0
