from utils.io import load_experimental_json, process_sample_records, save_processed, load_experimental_csv
from utils.compare import compare_dataset, save_comparison_results
from kinetics import BrubachModel
from utils.parameters import REACTOR


def test_csv_roundtrip_and_compare(tmp_path):
    # Load raw JSON records and process to DataFrame
    records = load_experimental_json("data/raw/sample_experiments.json")
    df = process_sample_records(records)
    assert not df.empty

    # save processed CSV
    processed_csv = tmp_path / "data_processed.csv"
    save_processed(df, str(processed_csv))
    assert processed_csv.exists()

    # load via CSV loader
    records_csv = load_experimental_csv(str(processed_csv))
    assert len(records_csv) == len(records)

    # run compare and save results
    model = BrubachModel()
    comp = compare_dataset(records, model, REACTOR, z_points=41)
    assert not comp.empty

    out_csv = tmp_path / "comparison_results.csv"
    save_comparison_results(comp, str(out_csv))
    assert out_csv.exists()
    # basic content checks
    dfcomp = comp
    assert "dCO_out" in dfcomp.columns
    assert "sel_model" in dfcomp.columns
