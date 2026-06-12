"""Unit tests for the severity flag (receptor-LoF + agonist -> futile)."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.config import Config  # noqa: E402
from repurpose import steps  # noqa: E402


def _cfg(enabled=True, damp=0.7):
    return Config(
        raw={"severity": {"enabled": enabled, "damp_factor": damp}},
        root=Path("."),
    )


def _ranked():
    return pd.DataFrame([
        # AGONIST drug on INSR + 'INSR deficiency' -> severe LoF trap
        {"rank": 1, "drug_id": "INS_HUMAN", "efo_id": "E1",
         "disease_name": "hyperinsulinism due to INSR deficiency",
         "disease_gene_match": "INSR", "opportunity": 1.0},
        # AGONIST drug on PPARG + 'PPARG-related' -> NOT 'deficiency' -> kept
        {"rank": 2, "drug_id": "PIO", "efo_id": "E2",
         "disease_name": "PPARG-related familial partial lipodystrophy",
         "disease_gene_match": "PPARG", "opportunity": 1.0},
        # No disease_gene_match -> never flagged
        {"rank": 3, "drug_id": "OTHER", "efo_id": "E3",
         "disease_name": "some other disease",
         "disease_gene_match": "", "opportunity": 1.0},
        # INHIBITOR drug on a 'deficiency' disease -> NOT flagged (inhibitor
        # of a broken receptor doesn't help either but that's a separate
        # concern; severity flag is specifically about agonist + LoF)
        {"rank": 4, "drug_id": "AB_X", "efo_id": "E4",
         "disease_name": "X deficiency syndrome",
         "disease_gene_match": "X", "opportunity": 1.0},
    ])


def _drug_targets():
    return pd.DataFrame([
        {"drug_id": "INS_HUMAN", "target_symbol": "INSR",  "action_type": "AGONIST",  "is_direct": True},
        {"drug_id": "PIO",       "target_symbol": "PPARG", "action_type": "AGONIST",  "is_direct": True},
        {"drug_id": "AB_X",      "target_symbol": "X",     "action_type": "INHIBITOR", "is_direct": True},
        {"drug_id": "OTHER",     "target_symbol": "Z",     "action_type": "AGONIST",  "is_direct": True},
    ])


def test_severity_flags_loF_agonist_only():
    out = steps.flag_severity_concern(_ranked(), _drug_targets(), _cfg())
    flagged = out[out.severity_concern == "severe_loF_agonist"]
    assert set(flagged["drug_id"]) == {"INS_HUMAN"}, \
        "must flag INSULIN -> INSR-deficiency; must NOT flag pioglitazone or X inhibitor"


def test_severity_damps_and_reranks():
    out = steps.flag_severity_concern(_ranked(), _drug_targets(), _cfg(damp=0.7))
    ins = out[out.drug_id == "INS_HUMAN"].iloc[0]
    # Damped 1.0 -> 0.3 (1 - 0.7)
    assert abs(ins["opportunity"] - 0.3) < 1e-4
    # And it should have fallen below the un-damped rows
    assert ins["rank"] > 1


def test_severity_inspection_only_when_damp_zero():
    """damp_factor=0 leaves opportunity unchanged but still populates the flag."""
    out = steps.flag_severity_concern(_ranked(), _drug_targets(), _cfg(damp=0.0))
    assert (out["severity_concern"] == "severe_loF_agonist").sum() == 1
    # Opportunities unchanged
    assert (out["opportunity"] == _ranked()["opportunity"].values).all()


def test_disabled_passes_through():
    out = steps.flag_severity_concern(_ranked(), _drug_targets(),
                                      _cfg(enabled=False))
    assert (out["severity_concern"] == "").all()


def test_partial_loF_disease_not_flagged():
    """'PPARG-related familial partial lipodystrophy' lacks the 'deficiency'
    keyword and must NOT be flagged, even though PPARG is the disease gene
    and TZD is an agonist."""
    out = steps.flag_severity_concern(_ranked(), _drug_targets(), _cfg())
    pio = out[out.drug_id == "PIO"].iloc[0]
    assert pio["severity_concern"] == ""
    # Opportunity untouched
    assert abs(pio["opportunity"] - 1.0) < 1e-4


def test_missing_columns_returns_unchanged():
    """If disease_gene_match isn't in the input, return the frame untouched
    (just with an empty severity_concern column)."""
    no_match = pd.DataFrame([{"drug_id": "X", "efo_id": "E1",
                              "disease_name": "X deficiency", "opportunity": 1.0}])
    out = steps.flag_severity_concern(no_match, _drug_targets(), _cfg())
    assert (out["severity_concern"] == "").all()


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"OK {name}")
