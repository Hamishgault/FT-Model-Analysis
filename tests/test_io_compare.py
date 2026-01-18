from utils.io import load_experimental_json, process_sample_records
from utils.compare import compare_dataset
from kinetics import BrubachModel
from utils.parameters import REACTOR


def test_process_and_compare_sample(tmp_path):
    records = load_experimental_json("data/raw/sample_experiments.json")
    df = process_sample_records(records)
    assert not df.empty

    # run a quick compare using the Brubach model (short z grid to be fast)
    model = BrubachModel()
    comp = compare_dataset(records, model, REACTOR, z_points=41)
    assert not comp.empty
    assert "dCO_out" in comp.columns
