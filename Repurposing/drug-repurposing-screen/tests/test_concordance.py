"""Unit tests for the disease-gene-concordance flag.

Flags (drug, disease) pairs where one of the drug's direct targets is
named in the disease itself -- the 'drug targets the broken receptor'
class. Inspection-only; no automatic score change.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose import steps  # noqa: E402


def _gene_info(symbols):
    return pd.DataFrame({"symbol": symbols, "gene_name": [""] * len(symbols)})


def test_extracts_direct_symbol():
    hyp = pd.DataFrame([{
        "drug_id": "D1", "efo_id": "E1",
        "disease_name": "hyperinsulinism due to INSR deficiency",
    }])
    dt = pd.DataFrame([{"drug_id": "D1", "target_symbol": "INSR",
                        "action_type": "AGONIST", "is_direct": True}])
    out = steps.flag_disease_gene_concordance(hyp, dt, _gene_info(["INSR"]))
    assert out.iloc[0]["disease_gene_match"] == "INSR"


def test_handles_receptor_suffix():
    """'X receptor' should derive XR (TSH receptor -> TSHR)."""
    hyp = pd.DataFrame([{
        "drug_id": "D1", "efo_id": "E1",
        "disease_name": "hypothyroidism due to TSH receptor mutations",
    }])
    dt = pd.DataFrame([{"drug_id": "D1", "target_symbol": "TSHR",
                        "action_type": "AGONIST", "is_direct": True}])
    out = steps.flag_disease_gene_concordance(hyp, dt, _gene_info(["TSHR"]))
    assert out.iloc[0]["disease_gene_match"] == "TSHR"


def test_partial_lof_legit_case_also_flags():
    """TZDs -> PPARG-related FPLD3 -- legit clinical strategy, but flag still fires
    (the flag is for human review, not auto-damping)."""
    hyp = pd.DataFrame([{
        "drug_id": "PIO", "efo_id": "E1",
        "disease_name": "PPARG-related familial partial lipodystrophy",
    }])
    dt = pd.DataFrame([{"drug_id": "PIO", "target_symbol": "PPARG",
                        "action_type": "AGONIST", "is_direct": True}])
    out = steps.flag_disease_gene_concordance(hyp, dt, _gene_info(["PPARG"]))
    assert out.iloc[0]["disease_gene_match"] == "PPARG"


def test_no_match_when_disease_mentions_other_gene():
    """Drug's targets aren't in the disease name -> empty match."""
    hyp = pd.DataFrame([{
        "drug_id": "D1", "efo_id": "E1",
        "disease_name": "hyperinsulinism due to INSR deficiency",
    }])
    dt = pd.DataFrame([{"drug_id": "D1", "target_symbol": "EGFR",
                        "action_type": "INHIBITOR", "is_direct": True}])
    out = steps.flag_disease_gene_concordance(hyp, dt, _gene_info(["INSR", "EGFR"]))
    assert out.iloc[0]["disease_gene_match"] == ""


def test_ignores_uppercase_non_gene_tokens():
    """'US' or 'EU' aren't genes -- must not be extracted as matches."""
    hyp = pd.DataFrame([{
        "drug_id": "D1", "efo_id": "E1",
        "disease_name": "US-approved indication for hypertension",
    }])
    dt = pd.DataFrame([{"drug_id": "D1", "target_symbol": "ADRA1A",
                        "action_type": "INHIBITOR", "is_direct": True}])
    out = steps.flag_disease_gene_concordance(hyp, dt, _gene_info(["ADRA1A"]))
    assert out.iloc[0]["disease_gene_match"] == ""


def test_excludes_string_neighbours():
    """Only direct targets are considered."""
    hyp = pd.DataFrame([{
        "drug_id": "D1", "efo_id": "E1",
        "disease_name": "INSR deficiency syndrome",
    }])
    dt = pd.DataFrame([
        {"drug_id": "D1", "target_symbol": "EGFR", "action_type": "INHIBITOR", "is_direct": True},
        # INSR shows up only as a STRING neighbour -- should NOT be reported
        {"drug_id": "D1", "target_symbol": "INSR", "action_type": "",          "is_direct": False},
    ])
    out = steps.flag_disease_gene_concordance(hyp, dt, _gene_info(["INSR", "EGFR"]))
    assert out.iloc[0]["disease_gene_match"] == ""


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"OK {name}")
