from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_orientation_table():
    df = pd.read_csv(ROOT / "inputs/gb_orientations.csv")
    assert len(df) == 1291
    assert df["id"].is_unique
    assert set(df["gb_type"].dropna()) <= {"tilt", "asymmetric_tilt", "twist", "mixed", "unclassified"}
