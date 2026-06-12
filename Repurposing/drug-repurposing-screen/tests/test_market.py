"""Unit tests for MarketSizeClient and add_market_size step."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.config import Config  # noqa: E402
from repurpose.sources.market import MarketSizeClient  # noqa: E402
from repurpose import steps  # noqa: E402


CURATED_HEADER = "efo_id,disease_name,us_patients,us_prevalence_per_100k,source,as_of\n"


def _csv(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "prev.csv"
    p.write_text(CURATED_HEADER + body)
    return p


def _cfg(top_n=10, orphan_threshold=200_000):
    return Config(
        raw={"market": {"top_n": top_n,
                        "rare_disease_us_threshold": orphan_threshold}},
        root=Path("."),
    )


def _ranked():
    return pd.DataFrame([
        {"rank": 1, "drug_id": "D1", "efo_id": "EFO_BC",
         "disease_name": "Breast cancer", "opportunity": 0.9, "lead_target": "T1"},
        {"rank": 2, "drug_id": "D2", "efo_id": "EFO_IPF",
         "disease_name": "IPF", "opportunity": 0.8, "lead_target": "T2"},
        {"rank": 3, "drug_id": "D3", "efo_id": "EFO_UNKNOWN",
         "disease_name": "Made up", "opportunity": 0.7, "lead_target": "T3"},
    ])


# ----------------------------------------------------------------------
# Client behaviour
# ----------------------------------------------------------------------
def test_lookup_basic(tmp_path):
    csv = _csv(tmp_path,
               "EFO_BC,breast cancer,4100000,1235,SEER,2023\n"
               "EFO_IPF,IPF,100000,30,PFF,2020\n")
    client = MarketSizeClient(curated_csv=csv)
    out = client.lookup(["EFO_BC", "EFO_IPF", "EFO_UNKNOWN"])
    assert len(out) == 3
    bc = out[out.efo_id == "EFO_BC"].iloc[0]
    ipf = out[out.efo_id == "EFO_IPF"].iloc[0]
    unk = out[out.efo_id == "EFO_UNKNOWN"].iloc[0]
    assert int(bc["us_patients"]) == 4_100_000
    assert bc["market_source"] == "SEER"
    assert int(ipf["us_patients"]) == 100_000
    # Unknown disease comes back with NaN/empty
    assert pd.isna(unk["us_patients"])
    assert unk["market_source"] == ""


def test_missing_csv_returns_empty_data(tmp_path):
    """Pipeline should still run if the curated CSV doesn't exist yet."""
    client = MarketSizeClient(curated_csv=tmp_path / "does_not_exist.csv")
    out = client.lookup(["EFO_BC"])
    assert pd.isna(out.iloc[0]["us_patients"])
    assert out.iloc[0]["market_source"] == ""


def test_missing_required_columns_raises(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("efo_id,us_patients\nEFO_BC,100\n")
    try:
        MarketSizeClient(curated_csv=bad)
    except KeyError as e:
        assert "us_prevalence_per_100k" in str(e) or "source" in str(e) or "as_of" in str(e)
    else:
        raise AssertionError("expected KeyError on malformed CSV")


def test_orphanet_fills_gaps_without_overriding_curated(tmp_path):
    """Curated takes precedence; Orphanet adds rows the curated CSV doesn't have."""
    cur = _csv(tmp_path, "EFO_BC,breast cancer,4100000,1235,SEER,2023\n")
    orph = tmp_path / "orph.csv"
    orph.write_text(
        CURATED_HEADER +
        "EFO_BC,breast cancer,9999,9.99,WRONG,2099\n"          # collision with curated
        "Orphanet_42,rare syndrome,166,0.05,Orphanet,2025\n"   # new rare disease
    )
    client = MarketSizeClient(curated_csv=cur, orphanet_csv=orph)
    out = client.lookup(["EFO_BC", "Orphanet_42", "EFO_UNKNOWN"])
    bc = out[out.efo_id == "EFO_BC"].iloc[0]
    # Curated wins on collision
    assert int(bc["us_patients"]) == 4_100_000
    assert bc["market_source"] == "SEER"
    rare = out[out.efo_id == "Orphanet_42"].iloc[0]
    assert int(rare["us_patients"]) == 166
    assert rare["market_source"] == "Orphanet"
    assert pd.isna(out[out.efo_id == "EFO_UNKNOWN"].iloc[0]["us_patients"])


def test_orphanet_only_works_without_curated(tmp_path):
    """Pure Orphanet backend (curated missing) still serves lookups."""
    orph = tmp_path / "orph.csv"
    orph.write_text(
        CURATED_HEADER +
        "Orphanet_42,rare syndrome,166,0.05,Orphanet,2025\n"
    )
    client = MarketSizeClient(curated_csv=None, orphanet_csv=orph)
    out = client.lookup(["Orphanet_42"])
    assert int(out.iloc[0]["us_patients"]) == 166


def test_duplicate_efo_id_first_wins(tmp_path):
    csv = _csv(tmp_path,
               "EFO_BC,breast cancer,4100000,1235,SEER,2023\n"
               "EFO_BC,breast cancer,9999999,9999,WRONG,2099\n")
    client = MarketSizeClient(curated_csv=csv)
    out = client.lookup(["EFO_BC"])
    assert int(out.iloc[0]["us_patients"]) == 4_100_000
    assert out.iloc[0]["market_source"] == "SEER"


# ----------------------------------------------------------------------
# Step integration -- add_market_size
# ----------------------------------------------------------------------
def test_step_attaches_market_columns(tmp_path):
    csv = _csv(tmp_path,
               "EFO_BC,breast cancer,4100000,1235,SEER,2023\n"
               "EFO_IPF,IPF,100000,30,PFF,2020\n")
    client = MarketSizeClient(curated_csv=csv)
    out = steps.add_market_size(_ranked(), client, _cfg())
    bc = out[out.drug_id == "D1"].iloc[0]
    ipf = out[out.drug_id == "D2"].iloc[0]
    unk = out[out.drug_id == "D3"].iloc[0]
    assert int(bc["us_patients"]) == 4_100_000
    assert int(ipf["us_patients"]) == 100_000
    # Orphan flag: IPF below 200k -> orphan; BC above -> not orphan; unknown -> NA
    assert ipf["is_orphan"] is True or ipf["is_orphan"] == True  # boolean dtype
    assert bc["is_orphan"] is False or bc["is_orphan"] == False
    assert pd.isna(unk["is_orphan"])


def test_step_does_not_change_ranking(tmp_path):
    """Market pass is inspection-only: never modifies opportunity or rank."""
    csv = _csv(tmp_path,
               "EFO_BC,breast cancer,4100000,1235,SEER,2023\n"
               "EFO_IPF,IPF,100000,30,PFF,2020\n")
    client = MarketSizeClient(curated_csv=csv)
    original = _ranked()
    out = steps.add_market_size(original, client, _cfg())
    assert out["rank"].tolist() == original["rank"].tolist()
    assert out["opportunity"].tolist() == original["opportunity"].tolist()
    assert out["drug_id"].tolist() == original["drug_id"].tolist()


def test_step_respects_top_n(tmp_path):
    """Rows outside top_n still get market columns (as NA), but no lookup happens."""
    csv = _csv(tmp_path,
               "EFO_BC,breast cancer,4100000,1235,SEER,2023\n"
               "EFO_IPF,IPF,100000,30,PFF,2020\n")
    client = MarketSizeClient(curated_csv=csv)
    # top_n=1 -> only the top hypothesis gets enriched
    out = steps.add_market_size(_ranked(), client, _cfg(top_n=1))
    # All output rows still have the columns
    for c in ["us_patients", "us_prevalence_per_100k", "market_source", "as_of", "is_orphan"]:
        assert c in out.columns
    # But only the top one (EFO_BC) was looked up
    bc = out[out.drug_id == "D1"].iloc[0]
    ipf = out[out.drug_id == "D2"].iloc[0]
    assert int(bc["us_patients"]) == 4_100_000
    assert pd.isna(ipf["us_patients"])


def test_step_handles_empty_ranked(tmp_path):
    csv = _csv(tmp_path, "EFO_BC,breast cancer,4100000,1235,SEER,2023\n")
    client = MarketSizeClient(curated_csv=csv)
    empty = pd.DataFrame(columns=["rank", "drug_id", "efo_id", "disease_name",
                                  "opportunity", "lead_target"])
    out = steps.add_market_size(empty, client, _cfg())
    assert len(out) == 0
    for c in ["us_patients", "us_prevalence_per_100k", "market_source", "as_of", "is_orphan"]:
        assert c in out.columns


def test_orphan_threshold_is_configurable(tmp_path):
    """Custom orphan threshold flips classification."""
    csv = _csv(tmp_path,
               "EFO_BC,breast cancer,4100000,1235,SEER,2023\n"
               "EFO_IPF,IPF,100000,30,PFF,2020\n")
    client = MarketSizeClient(curated_csv=csv)
    # Push the threshold above breast cancer's 4.1M -> everything is orphan now
    out = steps.add_market_size(_ranked(), client, _cfg(orphan_threshold=10_000_000))
    bc = out[out.drug_id == "D1"].iloc[0]
    assert bc["is_orphan"] is True or bc["is_orphan"] == True


if __name__ == "__main__":
    import tempfile
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            with tempfile.TemporaryDirectory() as d:
                fn(Path(d))
                print(f"OK {name}")
