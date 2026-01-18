"""I/O helpers for experimental datasets.

Functions here are small, focused utilities to load experimental data files
from `data/raw/` and produce standardized processed tables for model
comparison.
"""
from pathlib import Path
from typing import Dict, List

import pandas as pd


def load_experimental_json(path: str) -> List[Dict]:
    """Load an experimental dataset saved as JSON array of records."""
    p = Path(path)
    df = pd.read_json(p)
    return df.to_dict(orient="records")


def load_experimental_csv(path: str) -> List[Dict]:
    """Load an experimental dataset saved as CSV.

    Returns a list of records (dicts) compatible with the JSON loader output.
    """
    p = Path(path)
    df = pd.read_csv(p)
    return df.to_dict(orient="records")


def save_processed(df: pd.DataFrame, outpath: str):
    p = Path(outpath)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False)


def process_sample_records(records: List[Dict]) -> pd.DataFrame:
    """Convert raw records into a processed DataFrame with computed groups.

    The processed DataFrame includes computed product group fractions:
    - CH4 (HC_C1)
    - C2_4 (HC_C2 + HC_C3 + HC_C4)
    - C5+ (HC_C5 + HC_C6_10 + HC_C11_20 + HC_C21_30)

    Also keeps metadata fields like T, H2_CO2, GHSV.
    """
    rows = []
    for r in records:
        ch4 = r.get("HC_C1", 0.0)
        c2_4 = r.get("HC_C2", 0.0) + r.get("HC_C3", 0.0) + r.get("HC_C4", 0.0)
        c5p = r.get("HC_C5", 0.0) + r.get("HC_C6_10", 0.0) + r.get("HC_C11_20", 0.0) + r.get("HC_C21_30", 0.0)
        total = ch4 + c2_4 + c5p
        if total > 0:
            sel_ch4 = ch4 / total
            sel_c2_4 = c2_4 / total
            sel_c5p = c5p / total
        else:
            sel_ch4 = sel_c2_4 = sel_c5p = 0.0

        rows.append({
            "id": r.get("id"),
            "T": r.get("T"),
            "H2_CO2": r.get("H2_CO2"),
            "GHSV": r.get("GHSV"),
            "X_CO2": r.get("X_CO2"),
            "X_H2": r.get("X_H2"),
            "CO_out": r.get("CO_out"),
            "sel_CH4": sel_ch4,
            "sel_C2_4": sel_c2_4,
            "sel_C5_plus": sel_c5p,
        })

    return pd.DataFrame(rows)
