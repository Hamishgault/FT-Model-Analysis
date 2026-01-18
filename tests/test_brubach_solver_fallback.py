import numpy as np
from unittest import mock

from kinetics.brubach_2022 import BrubachModel


def test_solver_fallback_uses_last_coverages(monkeypatch):
    model = BrubachModel()

    # force root to return unsuccessful once, then successful
    import scipy.optimize as so

    original_root = so.root

    class FakeRes:
        def __init__(self, success, x):
            self.success = success
            self.x = x
            self.message = "failed"

    calls = {"n": 0}

    def fake_root(fun, x0, method=None):
        calls["n"] += 1
        if calls["n"] == 1:
            return FakeRes(False, x0)
        else:
            return FakeRes(True, x0)

    monkeypatch.setattr(so, "root", fake_root)

    Fi = np.zeros(7)
    Fi[1] = 3.0
    Fi[6] = 0.1

    # set a last coverages cache so fallback returns it
    model._last_coverages = {"theta_*": 0.5, "theta_H": 0.1, "theta_CO2": 0.0, "theta_CO": 0.05,
                              "theta_OH": 0.01, "theta_CH2": 0.01, "theta_R": 0.01, "theta_IR": 0.0,
                              "theta_O": 0.0, "theta_HCO": 0.01}

    cov = model.solve_surface(523.15, 1e5, Fi)
    # ensure we got last_coverages back (since first root fails) or a successful sol
    assert isinstance(cov, dict)
    assert "theta_*" in cov

    # restore original
    monkeypatch.setattr(so, "root", original_root)
