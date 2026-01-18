"""Reduced implementation of Brübach et al. (2022) CO2-FTS kinetic model.

This is a simplified, implementation-ready mechanistic model capturing the
main pathways (H2/CO2 adsorption equilibria, CO adsorption kinetics, RWGS,
H-assisted CO dissociation to CH2*, chain initiation/growth, and lumped
termination to CH4, C2-4, C5+). Surface coverages are solved using a
steady-state root-finding approach with a reduced set of coverages.

This model intentionally lumps chain species (R_n and IR_n) and termination
channels for initial reproducible behavior. It uses the parameter set
provided in the specification and exposes a per-model `nu` mapping so the
`PFR` can use model.nu directly.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
from scipy.optimize import root

from .base_model import KineticModel
from utils.species import SPECIES_IDX

R_gas = 8.31446261815324  # J/mol/K


@dataclass
class BrubachParams:
    K1: float = 5.82e-1  # bar^-1 (H2 adsorption equilibrium)
    K2: float = 1.37     # bar^-1 (CO2 adsorption equilibrium)
    k3p: float = 1.43e3  # CO adsorption forward, units consistent with Fi
    k3m: float = 9.21e3  # CO desorption
    k4p: float = 3.49e2  # H2O adsorption
    k4m: float = 1.45e1  # H2O desorption
    k5: float = 2.21e1   # RWGS effective (forward)
    k5m: float = 1e-2    # RWGS reverse (small default)
    k6: float = 4.10e1   # CO hydrogenation RDS (k6b)
    K6a: float = 1.0     # equilibrium constant for CO+H <-> HCO*
    k7: float = 1.64e-1  # initiation
    k8: float = 1.37e3   # chain growth
    k9a: float = 2.59e4  # methane termination (R1)
    k9b: float = 1.97e3  # alkane termination (R_n)
    k10a: float = 9.57e2
    k10b: float = 1.79e3
    k11a: float = 2.18e2
    k11b: float = 2.49e3
    k12: float = 1.0e2   # branched (iso) chain growth
    k13: float = 7.67e2  # iso-termination
    Gamma10: float = 0.452e3  # J/mol
    Gamma11b: float = 2.89e3  # J/mol
    phi: float = 1.0           # OH exponent in termination
    n_max: int = 20           # maximum chain length to resolve
    # catalyst loading (g catalyst per m^3 reactor volume) — used to convert
    # per-gram rates (mol g^-1 h^-1) to volumetric rates (mol m^-3 s^-1)
    # default catalyst loading set to 1e6 g/m^3 (1000 g per L) to represent
    # a typical packed-bed catalyst density on a reactor volume basis; adjust via
    # MODEL_PARAMS if needed
    cat_loading: float = 1e6  # g/m^3
    # activation energies are not used explicitly here but present for later
    EA4p: float = 151e3
    EA4m: float = 199e3
    EA10: float = 5.97e3


class BrubachModel(KineticModel):
    """Reduced Brubach 2022 kinetic model.

    State vector for surface coverages (algebraic unknowns):
        theta_* (free), theta_H, theta_CO2, theta_CO, theta_OH, theta_CH2, theta_R, theta_IR

    Reactions (returned rates) are ordered as:
        [r_CO_ads_des, r_RWGS, r_CO_to_CH2, r_initiation, r_chain_growth, r_term_CH4, r_term_C2_4, r_term_C5plus]

    The model supplies `nu` mapping these reactions to gas species
    [CO, H2, CH4, C2_4, C5plus, H2O, CO2].
    """

    def __init__(self, params: Optional[Dict] = None):
        self.params = BrubachParams(**(params or {})) if not isinstance(params, BrubachParams) else params
        self.n_rxns = 8

        # per-model stoichiometry mapping (rows: species order from utils.species SPECIES)
        # columns correspond to reactions returned by rate()
        # species order: CO, H2, CH4, C2_4, C5plus, H2O, CO2
        # r0: CO adsorption/desorption: CO + * <-> CO*
        # r1: RWGS: CO2 + H* -> CO + OH
        # r2: CO -> CH2 (consumes CO and H2)
        # r3: initiation (consumes CH2 and H)
        # r4: chain growth (net CH2 consumption)
        # r5: termination to CH4 (releases H2O)
        # r6: termination to C2_4
        # r7: termination to C5+
        self.nu = np.array([
            # r0  r1  r2  r3  r4  r5  r6  r7
            [-1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # CO
            [0.0, -1.0, -3.0, -1.0, 0.0, -3.0, -5.0, -12.0],  # H2 (approx H consumption equivalents)
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],   # CH4
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],   # C2_4
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],   # C5plus
            [0.0, 1.0, 1.0, 0.0, -1.0, 1.0, 2.0, 5.0],    # H2O (approx production; chain growth consumes CH2 so -1 in r4)
            [0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # CO2
        ])

    def solve_surface(self, T: float, P: float, Fi: np.ndarray):
        """Solve algebraic surface coverages (steady-state) for reduced set."""
        # Convert inlet molar flows Fi (mol/s) to mole fractions and partial pressures
        # Mole fractions
        F_total = max(Fi.sum(), 1e-12)
        y_CO = Fi[SPECIES_IDX["CO"]] / F_total
        y_H2 = Fi[SPECIES_IDX["H2"]] / F_total
        y_CO2 = Fi[SPECIES_IDX.get("CO2", -1)] / F_total if "CO2" in SPECIES_IDX else 0.0
        y_H2O = Fi[SPECIES_IDX["H2O"]] / F_total

        # Partial pressures in bar (parameters K1, K2 are in bar^-1)
        pCO_bar = y_CO * (P / 1e5)
        pH2_bar = y_H2 * (P / 1e5)
        pCO2_bar = y_CO2 * (P / 1e5)
        pH2O_bar = y_H2O * (P / 1e5)

        # Concentrations (mol/m^3) via ideal gas for volumetric rate expressions if needed
        c_total = P / (R_gas * T)
        cCO = y_CO * c_total
        cH2 = y_H2 * c_total
        cCO2 = y_CO2 * c_total
        cH2O = y_H2O * c_total

        # Use bar partial pressures for adsorption/equilibria and concentrations for volumetric terms
        p = {"CO": pCO_bar, "H2": pH2_bar, "CO2": pCO2_bar, "H2O": pH2O_bar}
        c = {"CO": cCO, "H2": cH2, "CO2": cCO2, "H2O": cH2O}
        prm = self.params

        # Unknowns: [theta_star, theta_H, theta_CO2, theta_CO, theta_OH, theta_CH2, theta_R, theta_IR, theta_O, theta_HCO]
        def resid(x):
            ths, thH, thCO2, thCO, thOH, thCH2, thR, thIR, thO, thHCO = x

            # adsorption equilibria (H2 and CO2 treated as quasi-equilibrium)
            eq_H = thH - ths * np.sqrt(prm.K1 * p["H2"])  # theta_H - theta_* sqrt(K1 pH2)
            eq_CO2 = thCO2 - prm.K2 * p["CO2"] * ths      # theta_CO2 - K2 pCO2 theta_*

            # CO adsorption/desorption kinetics
            r_co_ads = prm.k3p * p["CO"] * ths
            r_co_des = prm.k3m * thCO
            r0 = r_co_ads - r_co_des

            # H2O adsorption/desorption
            r4 = prm.k4p * p["H2O"] * ths**2 - prm.k4m * thOH * thH

            # RWGS (forward and small reverse)
            r_rwgs = prm.k5 * thCO2 * thH - prm.k5m * thCO * thOH * ths

            # HCO quasi-equilibrium: theta_HCO = K6a * theta_CO * theta_H
            eq_HCO = thHCO - prm.K6a * thCO * thH

            # RDS for CO hydrogenation (HCO + H -> ... -> CH2)
            r6b = prm.k6 * thHCO * thH

            # CO -> CH2 monomer formation is represented by r6b here
            r2 = r6b

            # initiation and chain growth/termination (lumped)
            r3 = prm.k7 * thCH2 * thH
            r4g = prm.k8 * thR * thCH2

            # lumped terminations to products (OH-promoted)
            r5 = prm.k9a * thR * thH * (thOH ** prm.phi)  # CH4
            r6 = prm.k9b * thR * thH * (thOH ** prm.phi)  # C2_4
            r7 = 0.1 * prm.k9b * thR * thH * (thOH ** prm.phi)  # C5+

            # branched growth and iso-termination (lumped)
            r11 = prm.k11a * thR * thCH2
            r12 = prm.k12 * thIR * thCH2
            r13 = prm.k13 * thIR * (thOH ** prm.phi)

            # Surface balances (net production = 0)
            # CO: produced by r0 (adsorption) and r_rwgs (reverse direction can consume/produce), consumed by HCO eq (implicit) and r2
            bal_CO = r0 + r_rwgs - r2

            # CH2: produced by r2, consumed by initiation and growth
            bal_CH2 = r2 - r3 - r4g

            # R: produced by initiation and growth, consumed by terminations and branching
            bal_R = r3 + r4g - (r5 + r6 + r7) - r11

            # OH: produced by RWGS and r4 (ads/des), consumed partially in terminations
            bal_OH = r_rwgs + r4 - (0.1 * (r5 + r6 + r7) + r13)

            # IR: branched pool balance
            bal_IR = r11 + r12 - r13

            # O* (surface oxygen) approximate balance: produced by CO2 dissociation (part of RWGS), consumed by OH formation (included in r_rwgs term)
            bal_O = -r_rwgs  # simplified

            # site balance
            bal_site = 1.0 - (ths + thH + thCO2 + thCO + thOH + thCH2 + thR + thIR + thO + thHCO)

            return np.array([eq_H, eq_CO2, bal_CO, bal_CH2, bal_R, bal_OH, bal_IR, bal_O, eq_HCO, bal_site])

        # Initial guess for the extended set
        x0 = np.array([0.5, 0.1, 0.0, 0.05, 0.01, 0.01, 0.01, 0.0, 0.0, 0.01])

        # Try a few strategies if root solver fails, and keep the last successful
        # coverages to use as initial guesses for the next call (continuation).
        try:
            sol = root(resid, x0, method="hybr")
        except Exception as e:
            sol = None
            logging = __import__("logging").getLogger(__name__)
            logging.debug("Initial root call raised: %s", e)

        # If initial attempt failed or returned unsuccessful, try robust fallbacks
        if (sol is None) or (not getattr(sol, "success", False)):
            logger = __import__("logging").getLogger(__name__)
            logger.debug("Surface solver did not converge on first attempt; trying fallbacks")

            # 1) try using last known coverages as initial guess
            if hasattr(self, "_last_coverages") and self._last_coverages is not None:
                x_last = np.array([
                    self._last_coverages.get("theta_*", 0.5),
                    self._last_coverages.get("theta_H", 0.1),
                    self._last_coverages.get("theta_CO2", 0.0),
                    self._last_coverages.get("theta_CO", 0.05),
                    self._last_coverages.get("theta_OH", 0.01),
                    self._last_coverages.get("theta_CH2", 0.01),
                    self._last_coverages.get("theta_R", 0.01),
                    self._last_coverages.get("theta_IR", 0.0),
                    self._last_coverages.get("theta_O", 0.0),
                    self._last_coverages.get("theta_HCO", 0.01),
                ])
                try:
                    sol = root(resid, x_last, method="hybr")
                except Exception:
                    sol = None

        if (sol is None) or (not getattr(sol, "success", False)):
            # 2) try least_squares with bounds (0,1)
            try:
                from scipy.optimize import least_squares

                lb = np.zeros_like(x0)
                ub = np.ones_like(x0)
                ls = least_squares(resid, x0, bounds=(lb, ub), xtol=1e-8, ftol=1e-8)
                if ls.success:
                    sol = ls
                    # provide compatibility: create object with x and success
                    class _S:
                        pass

                    s = _S()
                    s.x = ls.x
                    s.success = True
                    sol = s
            except Exception:
                sol = None

        if (sol is None) or (not getattr(sol, "success", False)):
            # final fallback: use last known coverages if available, or safe default
            logger = __import__("logging").getLogger(__name__)
            logger.warning("Surface solver failed after retries; using last known or default coverages")
            if hasattr(self, "_last_coverages") and self._last_coverages is not None:
                return dict(self._last_coverages)

            # safe default: mostly free sites
            coverages = {
                "theta_*": 0.98,
                "theta_H": 0.01,
                "theta_CO2": 0.0,
                "theta_CO": 0.0,
                "theta_OH": 0.0,
                # keep a tiny CH2 to allow chain initiation even on solver fallback
                "theta_CH2": 1e-4,
                # small R pool so termination/growth terms can be evaluated
                "theta_R": 1e-6,
                "theta_IR": 0.0,
                "theta_O": 0.0,
                "theta_HCO": 0.0,
            }
            # store fallback as last_coverages
            self._last_coverages = dict(coverages)
            return coverages

        # successful solution
        ths, thH, thCO2, thCO, thOH, thCH2, thR, thIR, thO, thHCO = sol.x
        coverages = {
            "theta_*": float(ths),
            "theta_H": float(thH),
            "theta_CO2": float(thCO2),
            "theta_CO": float(thCO),
            "theta_OH": float(thOH),
            "theta_CH2": float(thCH2),
            "theta_R": float(thR),
            "theta_IR": float(thIR),
            "theta_O": float(thO),
            "theta_HCO": float(thHCO),
        }

        # cache for continuation
        self._last_coverages = dict(coverages)
        return coverages

    def rate(self, T: float, P: float, Fi: np.ndarray) -> np.ndarray:
        # Solve surface coverages first
        cov = self.solve_surface(T, P, Fi)
        ths = cov["theta_*"]
        thH = cov["theta_H"]
        thCO2 = cov["theta_CO2"]
        thCO = cov["theta_CO"]
        thOH = cov["theta_OH"]
        thCH2 = cov["theta_CH2"]
        thR = cov["theta_R"]
        thO = cov.get("theta_O", 0.0)
        thHCO = cov.get("theta_HCO", 0.0)

        prm = self.params

        # rates consistent with solve_surface definitions (extended)
        # Recompute mole fractions / partial pressures for consistency
        F_total = max(Fi.sum(), 1e-12)
        y_CO = Fi[SPECIES_IDX["CO"]] / F_total
        y_H2 = Fi[SPECIES_IDX["H2"]] / F_total
        pCO_bar = y_CO * (P / 1e5)

        # adsorption/desorption uses partial pressure in bar
        r0 = prm.k3p * pCO_bar * ths - prm.k3m * thCO
        r_rwgs = prm.k5 * thCO2 * thH - prm.k5m * thCO * thOH * ths
        r2 = prm.k6 * thHCO * thH  # RDS CO hydrogenation via HCO* + H*
        r3 = prm.k7 * thCH2 * thH
        r4g = prm.k8 * thR * thCH2

        # termination channels
        r5 = prm.k9a * thR * thH * (thOH ** prm.phi)
        r6 = prm.k9b * thR * thH * (thOH ** prm.phi)
        r7 = 0.1 * prm.k9b * thR * thH * (thOH ** prm.phi)

        # Apply a temperature-dependent damping derived from Gamma10 to reduce
        # excessive short-chain termination (encourages chain growth at higher T)
        term_factor = np.exp(-prm.Gamma10 / (R_gas * T)) if prm.Gamma10 is not None else 1.0
        r5 *= term_factor
        r6 *= term_factor
        r7 *= term_factor

        # Interpret the model reaction rates as per-gram catalyst rates (mol g^-1 h^-1)
        # and convert to volumetric rates (mol m^-3 s^-1) using catalyst loading.
        # r [mol g^-1 h^-1] * cat_loading [g/m^3] / 3600 [s/h] -> mol m^-3 s^-1
        r = np.array([r0, r_rwgs, r2, r3, r4g, r5, r6, r7])
        r = r * (prm.cat_loading / 3600.0)

        return r
