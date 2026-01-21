"""Brubach et al. (2022) CO2-FTS kinetic model - Phase 1 (Corrected Reduced Implementation).

Reference: Brübach, L.; Hodonj, D.; Biffar, L.; Pfeifer, P. 
Detailed Kinetic Modeling of CO2-Based Fischer–Tropsch Synthesis. 
Catalysts 2022, 12(6), 630. https://doi.org/10.3390/catal12060630

This implementation uses the mechanisms and parameters from Table 4 of the paper:
- Direct CO2 dissociation for RWGS (Table 2)
- H-assisted CO dissociation for FTS (Table 3)
- All parameter values from paper's Table 4 (regression results)
- O* and HCO* treated implicitly via quasi-equilibrium (K5b, K6a)

Phase 1 (current): Corrected reduced model with proper mechanistic expressions
Phase 2 (future): Full chain-length tracking (θ_R_n, θ_IR_n individual unknowns)
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
    """Brubach et al. (2022) CO2-FTS kinetic parameters from Table 4.
    
    Reference conditions: T_ref = 300°C = 573.15 K, P = 10 bar.
    All rate constants in mol g^-1 h^-1 unless noted.
    All equilibrium constants dimensionless or per-bar as indicated.
    """
    K1: float = 5.82e-1  # bar^-1, H2 adsorption equilibrium (Table 4)
    K2: float = 1.37     # bar^-1, CO2 adsorption equilibrium (Table 4)
    k3p: float = 1.43e3  # mol g^-1 h^-1 bar^-1, CO adsorption (Table 4)
    k3m: float = 9.21e3  # mol g^-1 h^-1, CO desorption (Table 4)
    k4p: float = 3.49e2  # mol g^-1 h^-1 bar^-1, H2O adsorption (Table 4)
    k4m: float = 1.45e1  # mol g^-1 h^-1, H2O desorption (Table 4)
    k5: float = 2.21e1   # mol g^-1 h^-1, direct CO2 dissociation RWGS (Table 4)
    K6a: float = 4.10e1  # dimensionless, HCO* quasi-equilibrium constant (Table 4)
    k6: float = 1.64e-1  # mol g^-1 h^-1, H-assisted CO dissociation RDS (Table 4)
    k7: float = 9.57e2   # mol g^-1 h^-1, chain initiation (Table 4)
    k8: float = 1.79e3   # mol g^-1 h^-1, chain growth (Table 4)
    k9: float = 2.18e2   # mol g^-1 h^-1, n-alkane termination (Table 4)
    k10: float = 2.49e3  # mol g^-1 h^-1, 1-alkene termination (Table 4)
    k11: float = 7.67e2  # mol g^-1 h^-1, branching to iso-alkyl (Table 4)
    Gamma10: float = 4.52e2  # kJ/mol, chain-length dep for 1-alkene termination (Table 4)
    Gamma11: float = 2.89e3  # kJ/mol, chain-length dep for branching (Table 4)
    phi: float = 1.0           # OH exponent in termination damping (implicit in Table 1, r10)
    n_max: int = 20           # maximum chain length to resolve
    # catalyst loading (g catalyst per m^3 reactor volume) — used to convert
    # per-gram rates (mol g^-1 h^-1) to volumetric rates (mol m^-3 s^-1)
    cat_loading: float = 1e6  # g/m^3
    # activation energies (not yet implemented; for future temperature-dependent extension)
    EA4p: float = 151e3  # kJ/mol, H2O adsorption activation energy (Table 4)
    EA4m: float = 199e3  # kJ/mol, H2O desorption activation energy (Table 4)
    EA10: float = 5.97   # kJ/mol, 1-alkene termination activation energy (Table 4)


class BrubachModel(KineticModel):
    """Brubach et al. (2022) detailed kinetic model - Phase 1 (corrected reduced).

    Reference: Brübach et al., Catalysts 2022, 12(6), 630
    https://doi.org/10.3390/catal12060630

    Surface species (8 unknowns solved algebraically):
        θ_* (free site), θ_H, θ_CO2, θ_CO, θ_OH, θ_CH2, Σθ_R_n, Σθ_IR_n

    Key mechanism (from paper Tables 1-3):
    1. H2 + 2* ⇌ 2H*                (K1 equilibrium)
    2. CO2 + * ⇌ CO2*               (K2 equilibrium)
    3. CO + * ⇌ CO*                 (k3p forward, k3m reverse)
    4. H2O + 2* ⇌ OH* + H*          (k4p forward, k4m reverse)
    5. CO2* + H* → CO* + OH*        (RWGS, direct CO2 dissociation mechanism)
       Expression: r5 = k5 * θ_CO2 * θ_H / θ_OH
    6. CO* + 2H* → CH2* + OH*       (H-assisted CO dissociation, RDS via HCO*)
       Expression: r6 = k6 * K6a * θ_CO * θ_H^2 / θ_OH^phi
       (HCO* and O* are implicit via quasi-equilibrium per paper)
    7. CH2* + H* → R1*              (chain initiation)
    8. R_n* + CH2* → R_(n+1)*       (chain growth, constant α)
    9. R_n* + H* → P_n + *          (n-alkane termination)
    10. R_n* → Ol_n + H* + *        (1-alkene termination, β-hydride elimination)
    11. R_n* → IR_n*                (branching to iso-alkyl)
    12. IR_n* + CH2* → IR_(n+1)*    (iso-alkyl growth)
    13. IR_n* → I_n + H* + *        (iso-alkene termination)

    Reactions returned by rate():
        [r_CO_ads/des, r_RWGS, r_CO_hydrogenation, r_initiation, 
         r_chain_growth, r_term_alkane, r_term_alkene, r_branching]

    Stoichiometry nu maps these reactions to gas species [CO, H2, CH4, C2_4, C5plus, H2O, CO2].
    """

    def __init__(self, params: Optional[Dict] = None):
        self.params = BrubachParams(**(params or {})) if not isinstance(params, BrubachParams) else params
        self.n_rxns = 8

        # per-model stoichiometry (rows: species; columns: reactions)
        # species order: CO, H2, CH4, C2_4, C5plus, H2O, CO2
        self.nu = np.array([
            # r0  r1  r2  r3  r4  r5  r6  r7
            [-1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # CO
            [0.0, -1.0, -3.0, -1.0, 0.0, -3.0, -5.0, -12.0],  # H2
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],   # CH4
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],   # C2_4
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],   # C5plus
            [0.0, 1.0, 1.0, 0.0, -1.0, 1.0, 2.0, 5.0],  # H2O
            [0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # CO2
        ])

    def solve_surface(self, T: float, P: float, Fi: np.ndarray):
        """Solve algebraic surface coverages per Brubach et al. (2022).
        
        Uses direct CO2 dissociation mechanism (Table 2) with H-assisted CO dissociation (Table 3).
        O* and HCO* are implicit via quasi-equilibrium (K5b, K6a) — not tracked separately.
        
        Unknowns: [θ_*, θ_H, θ_CO2, θ_CO, θ_OH, θ_CH2, θ_R, θ_IR]
        Site balance: 1 = θ_* + θ_H + θ_CO2 + θ_CO + θ_OH + θ_CH2 + θ_R + θ_IR
        """
        # Convert inlet molar flows Fi (mol/s) to mole fractions and partial pressures
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

        p = {"CO": pCO_bar, "H2": pH2_bar, "CO2": pCO2_bar, "H2O": pH2O_bar}
        prm = self.params

        # Unknowns: [theta_*, theta_H, theta_CO2, theta_CO, theta_OH, theta_CH2, theta_R, theta_IR]
        def resid(x):
            ths, thH, thCO2, thCO, thOH, thCH2, thR, thIR = x

            # ===== Equilibrium expressions =====
            # Step 1: H2 adsorption (dissociative)
            # H2 + 2* ⇌ 2H*,  K1 = θ_H^2 / (θ_*^2 * pH2)
            eq_H = thH - ths * np.sqrt(prm.K1 * p["H2"])

            # Step 2: CO2 adsorption (associative)
            # CO2 + * ⇌ CO2*,  K2 = θ_CO2 / (θ_* * pCO2)
            eq_CO2 = thCO2 - prm.K2 * p["CO2"] * ths

            # ===== Kinetic rate expressions =====
            # Step 3: CO adsorption/desorption kinetics
            r_co_ads = prm.k3p * p["CO"] * ths
            r_co_des = prm.k3m * thCO
            r0 = r_co_ads - r_co_des

            # Step 4: H2O adsorption/desorption kinetics
            r4_ads = prm.k4p * p["H2O"] * ths**2
            r4_des = prm.k4m * thOH * thH
            r4 = r4_ads - r4_des

            # ===== RWGS Mechanism (Direct CO2 dissociation) =====
            # Step 5 (from Table 1): CO2* + H* → CO* + OH*
            # Quasi-equilibrium expression (O* implicit): r5 = k5 * θ_CO2 * θ_H / θ_OH
            thOH_safe = max(thOH, 1e-12)
            r5_rwgs = prm.k5 * thCO2 * thH / thOH_safe

            # ===== H-Assisted CO Dissociation (RDS for CO → CH2) =====
            # Steps 6a (quasi-eq) + 6b (RDS) from Table 1
            # CO* + 2H* → CH2* + OH*
            # Expression: r6 = k6 * K6a * θ_CO * θ_H^2 / θ_OH^phi
            thOH_phi = max(thOH ** prm.phi, 1e-12)
            r6_ch2 = prm.k6 * prm.K6a * thCO * (thH ** 2) / thOH_phi

            # ===== Chain Reactions =====
            # Step 7: Initiation (CH2* + H* → R1*)
            r7_init = prm.k7 * thCH2 * thH

            # Step 8: Growth (R_n* + CH2* → R_(n+1)*)
            r8_grow = prm.k8 * thR * thCH2

            # Step 9: n-alkane termination (R_n* + H* → P_n + *)
            r9_alkane = prm.k9 * thR * thH

            # Step 10: 1-alkene termination (R_n* → Ol_n + H* + *)
            # Expression from Table 1: k10 * exp(-Γ10*n/(RT)) * θ_R_n / θ_OH^phi
            # Simplified for lumped R: r10 = k10 * θ_R / θ_OH^phi
            r10_alkene = prm.k10 * thR / thOH_phi

            # Step 11: Branching (R_n* → IR_n*)
            # Expression: k11 * exp(-Γ11*n/(RT)) * θ_R_n * θ_CH2
            # Simplified: r11 = k11 * θ_R * θ_CH2
            r11_branch = prm.k11 * thR * thCH2

            # ISO-chain pseudo-reactions
            # Growth of ISO chains
            r12_iso_grow = 0.1 * prm.k11 * thIR * thCH2
            # Termination of ISO chains
            r13_iso_term = 0.1 * prm.k10 * thIR / thOH_phi

            # ===== Surface balance equations =====
            # Each balance: d(θ_k)/dt = Σ_j ν_k,j * r_j = 0 (steady-state)

            # Balance for CO*: produced by r0 and r5 (CO2 dissociation), consumed by r6 (hydrogenation)
            bal_CO = r0 + r5_rwgs - r6_ch2

            # Balance for CH2*: produced by r6, consumed by r7 (init) and r8 (growth) and r11 (branch)
            bal_CH2 = r6_ch2 - r7_init - r8_grow - r11_branch

            # Balance for R*: produced by r7 (init) and r8 (growth), consumed by r9, r10, r11
            bal_R = r7_init + r8_grow - r9_alkane - r10_alkene - r11_branch

            # Balance for OH*: produced by r5 (RWGS) and r4 (H2O ads)
            bal_OH = r5_rwgs + r4 - (0.01 * (r9_alkane + r10_alkene + r13_iso_term))

            # Balance for IR* (iso-alkyl chains)
            bal_IR = r11_branch + r12_iso_grow - r13_iso_term

            # Site balance (constrains θ_*)
            bal_site = 1.0 - (ths + thH + thCO2 + thCO + thOH + thCH2 + thR + thIR)

            return np.array([eq_H, eq_CO2, bal_CO, bal_CH2, bal_R, bal_OH, bal_IR, bal_site])

        # Initial guess
        x0 = np.array([0.5, 0.1, 0.05, 0.05, 0.05, 0.01, 0.01, 0.0])

        # Prefer bounded least_squares as primary solver
        sol = None

        # Use cached previous solution as initial guess for continuity
        if hasattr(self, "_last_coverages") and self._last_coverages is not None:
            x_init = np.array([
                self._last_coverages.get("theta_*", 0.5),
                self._last_coverages.get("theta_H", 0.1),
                self._last_coverages.get("theta_CO2", 0.05),
                self._last_coverages.get("theta_CO", 0.05),
                self._last_coverages.get("theta_OH", 0.05),
                self._last_coverages.get("theta_CH2", 0.01),
                self._last_coverages.get("theta_R", 0.01),
                self._last_coverages.get("theta_IR", 0.0),
            ])
        else:
            x_init = x0

        try:
            from scipy.optimize import least_squares

            lb = np.zeros_like(x0)
            ub = np.ones_like(x0)
            ls = least_squares(resid, x_init, bounds=(lb, ub), xtol=1e-8, ftol=1e-8, max_nfev=500)
            if ls.success:
                class _S:
                    pass
                s = _S()
                s.x = ls.x
                s.success = True
                sol = s
        except Exception:
            sol = None

        # Fallback to hybrid root solver
        if (sol is None) or (not getattr(sol, "success", False)):
            try:
                sol_root = root(resid, x_init, method="hybr")
                if getattr(sol_root, "success", False):
                    sol = sol_root
            except Exception:
                sol = None

        # Final fallback: use last known or safe defaults
        if (sol is None) or (not getattr(sol, "success", False)):
            logger = __import__("logging").getLogger(__name__)
            logger.warning("Surface solver failed; using last known or default coverages")
            if hasattr(self, "_last_coverages") and self._last_coverages is not None:
                return dict(self._last_coverages)

            # Safe default: mostly free sites with small nonzero CH2
            coverages = {
                "theta_*": 0.98,
                "theta_H": 0.01,
                "theta_CO2": 0.0,
                "theta_CO": 0.0,
                "theta_OH": 0.0,
                "theta_CH2": 1e-4,
                "theta_R": 1e-6,
                "theta_IR": 0.0,
            }
            self._last_coverages = dict(coverages)
            return coverages

        # Successful solution: extract and post-process
        ths, thH, thCO2, thCO, thOH, thCH2, thR, thIR = sol.x

        # Post-process: clamp negatives and normalize
        comps = {
            "theta_H": float(max(thH, 0.0)),
            "theta_CO2": float(max(thCO2, 0.0)),
            "theta_CO": float(max(thCO, 0.0)),
            "theta_OH": float(max(thOH, 0.0)),
            "theta_CH2": float(max(thCH2, 1e-10)),  # Keep nonzero for chain growth
            "theta_R": float(max(thR, 0.0)),
            "theta_IR": float(max(thIR, 0.0)),
        }

        # Normalize to ensure total occupancy ≤ 1
        sum_nonstar = sum(comps.values())
        if sum_nonstar >= 1.0:
            scale = 0.999 / sum_nonstar
            for k in comps:
                comps[k] *= scale
            ths = max(1e-6, 1.0 - sum(comps.values()))
        else:
            ths = max(1e-6, 1.0 - sum_nonstar)

        coverages = {"theta_*": float(ths)}
        coverages.update(comps)

        # Cache for continuation
        self._last_coverages = dict(coverages)
        return coverages

    def rate(self, T: float, P: float, Fi: np.ndarray) -> np.ndarray:
        """Compute reaction rates per Brubach et al. (2022) kinetic expressions.
        
        Returns rates for 8 reactions:
        [r0_CO_ads/des, r1_RWGS, r2_CO_hydrogenation, r3_initiation, 
         r4_growth, r5_alkane_term, r6_alkene_term, r7_branching]
        
        All rates converted from mol g^-1 h^-1 to mol m^-3 s^-1 using cat_loading.
        """
        # Solve surface coverages first
        cov = self.solve_surface(T, P, Fi)
        ths = cov["theta_*"]
        thH = cov["theta_H"]
        thCO2 = cov["theta_CO2"]
        thCO = cov["theta_CO"]
        thOH = cov["theta_OH"]
        thCH2 = cov["theta_CH2"]
        thR = cov["theta_R"]
        thIR = cov["theta_IR"]

        prm = self.params

        # Compute mole fractions / partial pressures
        F_total = max(Fi.sum(), 1e-12)
        y_CO = Fi[SPECIES_IDX["CO"]] / F_total
        y_H2 = Fi[SPECIES_IDX["H2"]] / F_total
        y_CO2 = Fi[SPECIES_IDX.get("CO2", -1)] / F_total if "CO2" in SPECIES_IDX else 0.0

        pCO_bar = y_CO * (P / 1e5)
        pH2_bar = y_H2 * (P / 1e5)
        pCO2_bar = y_CO2 * (P / 1e5)

        # ===== Compute reaction rates (mol g^-1 h^-1) =====

        # r0: CO adsorption/desorption (Step 3)
        r0 = prm.k3p * pCO_bar * ths - prm.k3m * thCO

        # r1: RWGS - direct CO2 dissociation (Step 5, Table 1)
        # Expression: r5 = k5 * θ_CO2 * θ_H / θ_OH
        thOH_safe = max(thOH, 1e-12)
        r1 = prm.k5 * thCO2 * thH / thOH_safe

        # r2: CO hydrogenation to CH2 (Steps 6a+6b, Table 1)
        # Expression: r6 = k6 * K6a * θ_CO * θ_H^2 / θ_OH^phi
        thOH_phi = max(thOH ** prm.phi, 1e-12)
        r2 = prm.k6 * prm.K6a * thCO * (thH ** 2) / thOH_phi

        # r3: Chain initiation (Step 7)
        # CH2* + H* → R1*
        r3 = prm.k7 * thCH2 * thH

        # r4: Chain growth (Step 8)
        # R_n* + CH2* → R_(n+1)*
        r4 = prm.k8 * thR * thCH2

        # r5: n-alkane termination (Step 9)
        # R_n* + H* → P_n + *
        r5 = prm.k9 * thR * thH

        # r6: 1-alkene termination (Step 10, β-hydride elimination)
        # Expression: r10 = k10 * θ_R / θ_OH^phi (simplified for lumped R)
        r6 = prm.k10 * thR / thOH_phi

        # r7: Branching to iso-alkyl (Step 11)
        # R_n* → IR_n*
        r7 = prm.k11 * thR * thCH2

        # Convert from mol g^-1 h^-1 to mol m^-3 s^-1
        # rate [mol g^-1 h^-1] * cat_loading [g/m^3] / 3600 [s/h] → mol m^-3 s^-1
        r = np.array([r0, r1, r2, r3, r4, r5, r6, r7])
        r = r * (prm.cat_loading / 3600.0)

        return r
