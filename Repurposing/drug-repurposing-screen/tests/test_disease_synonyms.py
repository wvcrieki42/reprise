"""Disease-synonym handling: adapter + loader + literature-pass integration.

OT exposes synonyms.hasExactSynonym as the high-precision variant; broad /
narrow / related synonyms are noisier and we exclude them.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.sources.loaders import load_disease_synonyms  # noqa: E402


def test_loader_missing_file_returns_empty(tmp_path):
    out = load_disease_synonyms(tmp_path / "does_not_exist.csv")
    assert out.empty
    assert list(out.columns) == ["efo_id", "disease_synonyms"]


def test_loader_reads_csv(tmp_path):
    p = tmp_path / "syn.csv"
    p.write_text(
        "efo_id,disease_synonyms\n"
        "EFO_0000305,carcinoma of breast;mammary carcinoma\n"
    )
    out = load_disease_synonyms(p)
    assert len(out) == 1
    assert out.iloc[0]["disease_synonyms"] == "carcinoma of breast;mammary carcinoma"


def test_loader_handles_missing_synonyms_column_gracefully(tmp_path):
    """Empty/NaN synonyms in a row should round-trip as empty string."""
    p = tmp_path / "syn.csv"
    p.write_text("efo_id,disease_synonyms\nEFO_X,\n")
    out = load_disease_synonyms(p)
    assert out.iloc[0]["disease_synonyms"] == ""


# ---------------------------------------------------------------- adapter
def _fake_diseases_parquet(tmp_path):
    """Write a tiny diseases parquet that mirrors OT's nested-synonym schema."""
    df = pd.DataFrame({
        "id": ["EFO_BC", "EFO_AD", "EFO_NONE"],
        "name": ["breast cancer", "Alzheimer disease", "made up syndrome"],
        "synonyms": [
            # exact + related + broad/narrow stub
            {"hasExactSynonym": np.array(["carcinoma of breast", "mammary carcinoma"]),
             "hasBroadSynonym": None,
             "hasNarrowSynonym": None,
             "hasRelatedSynonym": np.array(["BRCA"])},
            # The primary name MUST be excluded; case-insensitive dedup
            {"hasExactSynonym": np.array(["AD", "Alzheimer's disease",
                                          "alzheimer disease",     # case-only dup of primary
                                          "Alzheimer dementia"]),
             "hasBroadSynonym": None,
             "hasNarrowSynonym": None,
             "hasRelatedSynonym": None},
            # No exact synonyms -> row should be omitted
            {"hasExactSynonym": None,
             "hasBroadSynonym": np.array(["mystery"]),
             "hasNarrowSynonym": None,
             "hasRelatedSynonym": None},
        ],
    })
    parquet_dir = tmp_path / "diseases"
    parquet_dir.mkdir()
    df.to_parquet(parquet_dir / "part-00000.parquet")
    return parquet_dir


def test_adapter_keeps_exact_only_and_excludes_primary(tmp_path):
    from repurpose.sources import adapters
    diseases_dir = _fake_diseases_parquet(tmp_path)
    out = adapters.opentargets_disease_synonyms(str(diseases_dir))

    # Disease with only broad synonyms -> omitted
    assert "EFO_NONE" not in out["efo_id"].tolist()

    bc = out[out.efo_id == "EFO_BC"].iloc[0]
    assert bc["disease_synonyms"] == "carcinoma of breast;mammary carcinoma"
    # BRCA is a RELATED synonym -- must NOT appear (precision-focused).
    assert "BRCA" not in bc["disease_synonyms"]

    ad = out[out.efo_id == "EFO_AD"].iloc[0]
    syns = ad["disease_synonyms"].split(";")
    assert "AD" in syns
    assert "Alzheimer's disease" in syns
    assert "Alzheimer dementia" in syns
    # Case-insensitive duplicate of the primary disease name is excluded
    assert "alzheimer disease" not in [s.lower() for s in syns] or \
           syns.count("alzheimer disease") == 0


def test_adapter_respects_max_synonyms_cap(tmp_path):
    from repurpose.sources import adapters
    df = pd.DataFrame({
        "id": ["EFO_BIG"],
        "name": ["overloaded disease"],
        "synonyms": [{
            "hasExactSynonym": np.array([f"alias_{i}" for i in range(20)]),
            "hasBroadSynonym": None,
            "hasNarrowSynonym": None,
            "hasRelatedSynonym": None,
        }],
    })
    parquet_dir = tmp_path / "diseases"
    parquet_dir.mkdir()
    df.to_parquet(parquet_dir / "part.parquet")
    out = adapters.opentargets_disease_synonyms(str(parquet_dir), max_synonyms=3)
    syns = out.iloc[0]["disease_synonyms"].split(";")
    assert len(syns) == 3


if __name__ == "__main__":
    import tempfile
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            sig = fn.__code__.co_varnames[:fn.__code__.co_argcount]
            with tempfile.TemporaryDirectory() as d:
                if "tmp_path" in sig:
                    fn(Path(d))
                else:
                    fn()
                print(f"OK {name}")
