"""Unit tests for active-substance grouping.

Collapses (drug_id, efo_id) rows into (substance_chembl_id, efo_id)
keeping the highest-opportunity row per substance-disease as the
representative; preserves the variant names + count.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.config import Config  # noqa: E402
from repurpose import steps  # noqa: E402


def _cfg(enabled=True, collapse=True):
    return Config(
        raw={"substance_grouping": {"enabled": enabled, "collapse": collapse}},
        root=Path("."),
    )


def _ranked():
    """Three Pioglitazone variants + a different drug, all for the same disease."""
    return pd.DataFrame([
        {"rank": 1, "drug_id": "CHEMBL595",  "drug_name": "PIOGLITAZONE",
         "efo_id": "EFO_1", "opportunity": 1.05},
        {"rank": 2, "drug_id": "CHEMBL1715", "drug_name": "PIOGLITAZONE HYDROCHLORIDE",
         "efo_id": "EFO_1", "opportunity": 1.05},
        {"rank": 3, "drug_id": "CHEMBL3000", "drug_name": "ROSIGLITAZONE",
         "efo_id": "EFO_1", "opportunity": 1.05},
        # Same drug, different disease -- should stay separate
        {"rank": 4, "drug_id": "CHEMBL595",  "drug_name": "PIOGLITAZONE",
         "efo_id": "EFO_2", "opportunity": 0.8},
    ])


def _map():
    return pd.DataFrame([
        {"drug_id": "CHEMBL595",  "substance_chembl_id": "CHEMBL595",  "substance_name": "PIOGLITAZONE"},
        {"drug_id": "CHEMBL1715", "substance_chembl_id": "CHEMBL595",  "substance_name": "PIOGLITAZONE"},
        {"drug_id": "CHEMBL3000", "substance_chembl_id": "CHEMBL3000", "substance_name": "ROSIGLITAZONE"},
    ])


def test_collapse_groups_salt_with_parent():
    out = steps.collapse_to_substances(_ranked(), _map(), _cfg())
    # 3 variants on EFO_1 collapse to 2 substances (PIO + ROSI) + 1 PIO on EFO_2
    assert len(out) == 3
    pio_e1 = out[(out.substance_chembl_id == "CHEMBL595") & (out.efo_id == "EFO_1")].iloc[0]
    assert pio_e1["n_variants"] == 2
    assert "PIOGLITAZONE" in pio_e1["variant_names"]
    assert "PIOGLITAZONE HYDROCHLORIDE" in pio_e1["variant_names"]
    # Displayed name uses the canonical substance name, not the salt form
    assert pio_e1["drug_name"] == "PIOGLITAZONE"


def test_collapse_keeps_substances_for_other_diseases():
    """Same substance for a DIFFERENT disease stays as its own row."""
    out = steps.collapse_to_substances(_ranked(), _map(), _cfg())
    pio_e2 = out[(out.substance_chembl_id == "CHEMBL595") & (out.efo_id == "EFO_2")]
    assert len(pio_e2) == 1
    assert pio_e2.iloc[0]["n_variants"] == 1


def test_collapse_reranks_after_grouping():
    """Re-rank after collapse so the rank column is 1..N over the new row count."""
    out = steps.collapse_to_substances(_ranked(), _map(), _cfg())
    assert sorted(out["rank"].tolist()) == [1, 2, 3]


def test_collapse_can_be_disabled():
    """enabled=False -> pass through with no substance columns added."""
    out = steps.collapse_to_substances(_ranked(), _map(), _cfg(enabled=False))
    assert "substance_chembl_id" not in out.columns
    assert len(out) == 4  # no collapse


def test_columns_only_mode_skips_row_collapse():
    """enabled=True, collapse=False -> add substance columns but don't collapse rows."""
    out = steps.collapse_to_substances(_ranked(), _map(), _cfg(collapse=False))
    assert "substance_chembl_id" in out.columns
    assert len(out) == 4  # rows preserved
    # The two pioglitazone variants both have the canonical substance
    pio_rows = out[out.substance_chembl_id == "CHEMBL595"]
    assert len(pio_rows) == 3
    assert set(pio_rows["drug_name"]) == {"PIOGLITAZONE", "PIOGLITAZONE HYDROCHLORIDE"}


def test_empty_substance_map_passes_through():
    """When the substance map file is missing, every drug is its own substance."""
    empty = pd.DataFrame(columns=["drug_id", "substance_chembl_id", "substance_name"])
    out = steps.collapse_to_substances(_ranked(), empty, _cfg())
    assert len(out) == 4  # no collapse possible
    assert (out["substance_chembl_id"] == out["drug_id"]).all()


def test_highest_opportunity_wins_as_representative():
    """When variants have different opportunity scores, the highest wins."""
    ranked = pd.DataFrame([
        {"rank": 1, "drug_id": "CHEMBL595",  "drug_name": "PIOGLITAZONE",
         "efo_id": "EFO_1", "opportunity": 1.2},
        {"rank": 2, "drug_id": "CHEMBL1715", "drug_name": "PIOGLITAZONE HYDROCHLORIDE",
         "efo_id": "EFO_1", "opportunity": 0.9},
    ])
    out = steps.collapse_to_substances(ranked, _map(), _cfg())
    assert len(out) == 1
    assert out.iloc[0]["opportunity"] == 1.2


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"OK {name}")
