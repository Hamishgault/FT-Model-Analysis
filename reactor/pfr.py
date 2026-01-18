"""Generic 1D plug-flow reactor (PFR) solver.

The PFR integrates material balances:
    dF_i/dz = A * sum_j nu_ij * r_j(T, P, F)

where r_j are reaction rates [mol m^-3 s^-1].
"""
from typing import Optional
import numpy as np
from scipy.integrate import solve_ivp

from kinetics.base_model import KineticModel


class PFR:
    """Simple isothermal PFR solver.

    Parameters
    ----------
    A : float
        Cross-sectional area [m^2]
    L : float
        Reactor length [m]
    model : KineticModel
        Kinetic model instance implementing `rate(T,P,Fi)`
    nu : np.ndarray
        Stoichiometric coefficients matrix (n_species x n_rxns). Row i,
        column j is nu_ij.
    """

    def __init__(self, A: float, L: float, model: KineticModel, nu: Optional[np.ndarray] = None):
        self.A = A
        self.L = L
        self.model = model

        # Prefer the model's own nu (if provided). Otherwise use the nu passed
        # to the reactor. If neither is available, error out.
        if getattr(model, "nu", None) is not None:
            self.nu = np.asarray(model.nu)
        elif nu is not None:
            self.nu = np.asarray(nu)
        else:
            raise ValueError("No stoichiometric matrix `nu` provided. Either pass `nu` to PFR or provide `model.nu`.")

    def run(self, T: float, P: float, F0: np.ndarray, z_eval: Optional[np.ndarray] = None, **kwargs):
        """Integrate the PFR.

        Returns a `scipy.integrate.OdeResult` with `t` (z positions) and `y` (F_i).
        """
        F0 = np.asarray(F0)

        def odes(z, F):
            r = self.model.rate(T, P, F)
            r = np.atleast_1d(r)

            # If model returns a single overall rate but `nu` contains
            # multiple reaction columns, interpret the single rate as the
            # overall FT rate and distribute it across species according
            # to the column-sum of `nu` (i.e., treat it as one effective
            # lumped reaction). This keeps simple placeholder models
            # compatible with the fixed stoichiometric matrix.
            if r.size == 1 and self.nu.shape[1] != 1:
                nu_eff = self.nu.sum(axis=1)  # shape (n_species,)
                dFdz = self.A * (nu_eff * r[0])
                return dFdz

            if r.size != self.nu.shape[1]:
                raise ValueError(
                    f"Mismatch between number of rates (got {r.size}) and nu columns ({self.nu.shape[1]})."
                    " Provide a nu that matches the model or have the model return a rate vector matching nu."
                )

            dFdz = self.A * self.nu.dot(r)  # shape (n_species,)
            return dFdz

        z_span = (0.0, self.L)
        if z_eval is None:
            z_eval = np.linspace(0.0, self.L, 101)

        sol = solve_ivp(odes, z_span, F0, t_eval=z_eval, vectorized=False, **kwargs)
        return sol
