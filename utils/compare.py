"""Comparison helpers to match experimental runs with model outputs."""
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

from reactor.pfr import PFR
from utils.species import SPECIES_IDX


def make_inlet_from_experiment(rec: Dict, co2_basis: float = 1.0) -> np.ndarray:
    """Create inlet molar flows vector `Fi` from an experimental record.

    Convention: use CO2 as basis species (common in CO2-FTS experiments).
    `co2_basis` sets the inlet CO2 molar flow (default 1.0). H2 inlet is set
    via `H2_CO2` ratio in the record: H2 = H2_CO2 * co2_basis. CO inlet is set
    to zero (CO2 feed case).
    """
    n = len(SPECIES_IDX)
    Fi = np.zeros(n)
    Fi[SPECIES_IDX["CO2"]] = co2_basis
    Fi[SPECIES_IDX["H2"]] = float(rec.get("H2_CO2", 3.0)) * co2_basis
    Fi[SPECIES_IDX["CO"]] = 0.0
    return Fi


def compare_record_to_model(rec: Dict, model, reactor_params: Dict, z_points: int = 201) -> Dict:
    Fi = make_inlet_from_experiment(rec)
    pfr = PFR(reactor_params["A"], reactor_params["L"], model, nu=None)
    z_eval = np.linspace(0.0, reactor_params["L"], z_points)
    sol = pfr.run(T=rec.get("T", 523.15), P=1e5, F0=Fi, z_eval=z_eval)
    F_out = sol.y[:, -1]

    # model selectivity: use utils.species.compute_selectivity style (CH4, C2_4, C5+)
    iCH4 = SPECIES_IDX["CH4"]
    iC2_4 = SPECIES_IDX["C2_4"]
    iC5 = SPECIES_IDX["C5plus"]

    prod_CH4 = max(F_out[iCH4] - Fi[iCH4], 0.0)
    prod_C2_4 = max(F_out[iC2_4] - Fi[iC2_4], 0.0)
    prod_C5 = max(F_out[iC5] - Fi[iC5], 0.0)
    total = prod_CH4 + prod_C2_4 + prod_C5
    if total <= 0:
        sel_model = {"CH4": 0.0, "C2_4": 0.0, "C5+": 0.0}
    else:
        sel_model = {"CH4": prod_CH4 / total, "C2_4": prod_C2_4 / total, "C5+": prod_C5 / total}

    # experimental selectivity provided in processed data: sel_CH4, sel_C2_4, sel_C5_plus
    sel_exp = {"CH4": rec.get("sel_CH4", 0.0), "C2_4": rec.get("sel_C2_4", 0.0), "C5+": rec.get("sel_C5_plus", 0.0)}

    # compute simple errors
    errors = {
        "CO_out_model": float(F_out[SPECIES_IDX["CO"]]),
        "CO_out_exp": rec.get("CO_out"),
        "dCO_out": float(F_out[SPECIES_IDX["CO"]]) - float(rec.get("CO_out", 0.0)),
        "sel_model": sel_model,
        "sel_exp": sel_exp,
        "sel_error": {k: sel_model[k] - sel_exp.get(k, 0.0) for k in sel_model},
    }
    return errors


def compare_dataset(records: Iterable[Dict], model, reactor_params: Dict, z_points: int = 201) -> pd.DataFrame:
    rows = []
    for rec in records:
        r = compare_record_to_model(rec, model, reactor_params, z_points=z_points)
        rows.append({"id": rec.get("id"), "T": rec.get("T"), "H2_CO2": rec.get("H2_CO2"), **r})
    return pd.DataFrame(rows)


def save_comparison_results(df: pd.DataFrame, outpath: str):
    p = Path(outpath)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False)
