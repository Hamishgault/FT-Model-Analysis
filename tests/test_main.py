import matplotlib
matplotlib.use("Agg")  # non-interactive backend for tests

from pathlib import Path
import json

from main import run_and_report


def test_run_and_report_saves(tmp_path: Path):
    outdir = tmp_path / "out"
    # only Brubach is available in current configuration
    results, summaries = run_and_report(["Brubach"], save=True, outdir=outdir, z_points=51)

    assert "Brubach" in summaries

    conv = outdir / "conversion_co.png"
    sel = outdir / "selectivity.png"
    summ = outdir / "summaries.json"

    assert conv.exists()
    assert sel.exists()
    assert summ.exists()

    with open(summ) as fh:
        data = json.load(fh)
    assert "Brubach" in data
