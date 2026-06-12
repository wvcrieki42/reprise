"""Loader normalisation tests.

The repurposing pipeline keys ontology IDs in underscore form
(EFO_0001360, MONDO_0005148, HP_0000863). ChEMBL exports drug
indications with non-EFO sources in colon form (MONDO:0005148,
HP:0000863). Without normalisation, the novelty step silently fails
to subtract known indications for any drug whose ChEMBL indication is
coded against MONDO/HP -- approved drugs (SGLT2 inhibitors for T2D,
CFTR modulators for cystic fibrosis, etc.) end up flagged "novel".
This test locks the normalisation in.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.config import Config  # noqa: E402
from repurpose.sources.loaders import load_drug_indications  # noqa: E402
from repurpose import steps  # noqa: E402


def test_indications_normalise_colon_to_underscore(tmp_path):
    p = tmp_path / "indications.csv"
    p.write_text(
        "drug_id,efo_id,indication_name\n"
        "CHEMBL1,EFO_0000400,diabetes mellitus\n"        # already correct
        "CHEMBL2,MONDO:0005148,type 2 diabetes mellitus\n"  # colon -> underscore
        "CHEMBL3,HP:0000863,Central diabetes insipidus\n"   # colon -> underscore
        "CHEMBL4,DOID:0050700,liver cancer\n"               # colon -> underscore
    )
    df = load_drug_indications(p)
    ids = set(df["efo_id"])
    assert ids == {"EFO_0000400", "MONDO_0005148", "HP_0000863", "DOID_0050700"}
    # No colons survive anywhere in the column
    assert not df["efo_id"].str.contains(":", regex=False).any()


def test_indications_missing_indication_name_is_added(tmp_path):
    p = tmp_path / "indications.csv"
    p.write_text("drug_id,efo_id\nCHEMBL1,EFO_0000400\n")
    df = load_drug_indications(p)
    assert "indication_name" in df.columns
    assert df.iloc[0]["indication_name"] == ""


# ----------------------------------------------------------------------
# Target-class rollup -- closes ChEMBL indication-coverage gaps for
# formulations that share an active ingredient.
# ----------------------------------------------------------------------
def _cfg(enabled=True):
    return Config(raw={"novelty": {"expand_by_target_class": enabled}},
                  root=Path("."))


def test_target_class_rollup_pools_indications():
    """Two drugs with identical (target, action) sets share indications."""
    indications = pd.DataFrame([
        # INSULIN HUMAN has T2D + diabetes; the protamine variant has none
        {"drug_id": "INS_HUMAN",   "efo_id": "MONDO_0005148",
         "indication_name": "type 2 diabetes mellitus"},
        {"drug_id": "INS_HUMAN",   "efo_id": "EFO_0000400",
         "indication_name": "diabetes mellitus"},
    ])
    drug_targets = pd.DataFrame([
        {"drug_id": "INS_HUMAN",   "target_symbol": "INSR", "action_type": "AGONIST", "is_direct": True},
        {"drug_id": "INS_PROT",    "target_symbol": "INSR", "action_type": "AGONIST", "is_direct": True},
        # Different class: kinase inhibitor, must NOT pick up insulin's T2D
        {"drug_id": "OTHER_KIN",   "target_symbol": "EGFR", "action_type": "INHIBITOR", "is_direct": True},
    ])
    out = steps.expand_indications_by_target_class(indications, drug_targets, _cfg())
    # The protamine variant inherits BOTH of INSULIN HUMAN's indications
    prot = out[out.drug_id == "INS_PROT"]
    assert set(prot["efo_id"]) == {"MONDO_0005148", "EFO_0000400"}
    # The kinase inhibitor inherits nothing (different class)
    assert (out.drug_id == "OTHER_KIN").sum() == 0


def test_target_class_rollup_respects_action_type():
    """Same target, different action_type -> different class."""
    indications = pd.DataFrame([
        {"drug_id": "AGONIST_DRUG", "efo_id": "EFO_X", "indication_name": "x"},
    ])
    drug_targets = pd.DataFrame([
        {"drug_id": "AGONIST_DRUG", "target_symbol": "INSR", "action_type": "AGONIST", "is_direct": True},
        {"drug_id": "INHIB_DRUG",   "target_symbol": "INSR", "action_type": "INHIBITOR", "is_direct": True},
    ])
    out = steps.expand_indications_by_target_class(indications, drug_targets, _cfg())
    # INHIB_DRUG must NOT inherit AGONIST_DRUG's indication
    assert (out.drug_id == "INHIB_DRUG").sum() == 0


def test_target_class_rollup_excludes_string_neighbours():
    """STRING-expanded neighbours don't count toward the class key."""
    indications = pd.DataFrame([
        {"drug_id": "DRUG_A", "efo_id": "EFO_X", "indication_name": "x"},
    ])
    drug_targets = pd.DataFrame([
        {"drug_id": "DRUG_A", "target_symbol": "INSR", "action_type": "AGONIST", "is_direct": True},
        {"drug_id": "DRUG_A", "target_symbol": "IRS1", "action_type": "",        "is_direct": False},
        # DRUG_B has the SAME direct target but a different neighbour set;
        # rollup should still pool (neighbours don't change class identity)
        {"drug_id": "DRUG_B", "target_symbol": "INSR", "action_type": "AGONIST", "is_direct": True},
        {"drug_id": "DRUG_B", "target_symbol": "AKT1", "action_type": "",        "is_direct": False},
    ])
    out = steps.expand_indications_by_target_class(indications, drug_targets, _cfg())
    assert (out.drug_id == "DRUG_B").sum() == 1


def test_target_class_rollup_selective_inherits_from_broad():
    """A SELECTIVE drug (single target) inherits from a BROAD drug whose target
    set is a superset. Silodosin (ADRA1A-only) inherits doxazosin's (pan-ADRA1)
    hypertension indication so the screen doesn't surface 'use alpha-1 blocker
    for HTN' as a novel hypothesis."""
    indications = pd.DataFrame([
        {"drug_id": "DOX_BROAD", "efo_id": "EFO_HTN",
         "indication_name": "hypertension"},
    ])
    drug_targets = pd.DataFrame([
        # Selective: only ADRA1A
        {"drug_id": "SIL_SELECTIVE", "target_symbol": "ADRA1A",
         "action_type": "ANTAGONIST", "is_direct": True},
        # Broad: ADRA1A + ADRA1B + ADRA1D
        {"drug_id": "DOX_BROAD", "target_symbol": "ADRA1A",
         "action_type": "ANTAGONIST", "is_direct": True},
        {"drug_id": "DOX_BROAD", "target_symbol": "ADRA1B",
         "action_type": "ANTAGONIST", "is_direct": True},
        {"drug_id": "DOX_BROAD", "target_symbol": "ADRA1D",
         "action_type": "ANTAGONIST", "is_direct": True},
    ])
    out = steps.expand_indications_by_target_class(indications, drug_targets, _cfg())
    sil = out[out.drug_id == "SIL_SELECTIVE"]
    assert set(sil["efo_id"]) == {"EFO_HTN"}, \
        "selective drug must inherit broad drug's indication when its class is a subset"


def test_target_class_rollup_broad_does_not_inherit_from_selective():
    """The reverse: a BROADER drug does NOT inherit from a more selective one,
    because the selective drug's class is not a superset of the broad drug's."""
    indications = pd.DataFrame([
        {"drug_id": "SIL_SELECTIVE", "efo_id": "EFO_BPH",
         "indication_name": "benign prostatic hyperplasia"},
    ])
    drug_targets = pd.DataFrame([
        {"drug_id": "SIL_SELECTIVE", "target_symbol": "ADRA1A",
         "action_type": "ANTAGONIST", "is_direct": True},
        {"drug_id": "DOX_BROAD", "target_symbol": "ADRA1A",
         "action_type": "ANTAGONIST", "is_direct": True},
        {"drug_id": "DOX_BROAD", "target_symbol": "ADRA1B",
         "action_type": "ANTAGONIST", "is_direct": True},
    ])
    out = steps.expand_indications_by_target_class(indications, drug_targets, _cfg())
    # SIL has its own indication; DOX must NOT inherit BPH via the rollup
    assert (out.drug_id == "DOX_BROAD").sum() == 0


def test_target_class_rollup_unrelated_classes_dont_pool():
    """Two classes that aren't subset-related don't share indications."""
    indications = pd.DataFrame([
        {"drug_id": "DRUG_A", "efo_id": "EFO_X", "indication_name": "x"},
        {"drug_id": "DRUG_B", "efo_id": "EFO_Y", "indication_name": "y"},
    ])
    drug_targets = pd.DataFrame([
        # Different targets entirely
        {"drug_id": "DRUG_A", "target_symbol": "EGFR",
         "action_type": "INHIBITOR", "is_direct": True},
        {"drug_id": "DRUG_B", "target_symbol": "BRAF",
         "action_type": "INHIBITOR", "is_direct": True},
    ])
    out = steps.expand_indications_by_target_class(indications, drug_targets, _cfg())
    a = out[out.drug_id == "DRUG_A"]
    b = out[out.drug_id == "DRUG_B"]
    assert set(a["efo_id"]) == {"EFO_X"}
    assert set(b["efo_id"]) == {"EFO_Y"}


def test_target_class_rollup_can_be_disabled():
    """When disabled, indications pass through untouched."""
    indications = pd.DataFrame([
        {"drug_id": "DRUG_A", "efo_id": "EFO_X", "indication_name": "x"},
    ])
    drug_targets = pd.DataFrame([
        {"drug_id": "DRUG_A", "target_symbol": "INSR", "action_type": "AGONIST", "is_direct": True},
        {"drug_id": "DRUG_B", "target_symbol": "INSR", "action_type": "AGONIST", "is_direct": True},
    ])
    out = steps.expand_indications_by_target_class(indications, drug_targets, _cfg(enabled=False))
    assert out.equals(indications)


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
