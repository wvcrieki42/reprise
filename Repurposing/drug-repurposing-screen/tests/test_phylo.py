"""Unit tests for the phylogenetics (orthologous-gene model-organism) filter.

Boost-only design: present evidence multiplies opportunity by up to
(1 + boost_factor); absent evidence applies no penalty.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.config import Config  # noqa: E402
from repurpose import steps  # noqa: E402


def _cfg(boost=0.5):
    return Config(raw={"phylogenetics": {"boost_factor": boost}}, root=Path("."))


def _edges():
    """Synthetic drug -> target -> disease edges."""
    return pd.DataFrame([
        # D1 hits TGT_PHYLO_HIGH and TGT_PHYLO_LOW for E1 (novel disease)
        {"drug_id": "D1", "target_symbol": "TGT_PHYLO_HIGH", "efo_id": "E1"},
        {"drug_id": "D1", "target_symbol": "TGT_PHYLO_LOW",  "efo_id": "E1"},
        # D2 hits TGT_NO_EVIDENCE for E1
        {"drug_id": "D2", "target_symbol": "TGT_NO_EVIDENCE", "efo_id": "E1"},
        # D3 hits two targets for E2; only one has evidence
        {"drug_id": "D3", "target_symbol": "TGT_NO_EVIDENCE", "efo_id": "E2"},
        {"drug_id": "D3", "target_symbol": "TGT_PHYLO_MED",   "efo_id": "E2"},
    ])


def _phylo():
    return pd.DataFrame([
        {"target_symbol": "TGT_PHYLO_HIGH", "efo_id": "E1",
         "phylo_score": 1.0, "n_models": 5, "sources": "impc"},
        {"target_symbol": "TGT_PHYLO_LOW",  "efo_id": "E1",
         "phylo_score": 0.3, "n_models": 1, "sources": "impc"},
        {"target_symbol": "TGT_PHYLO_MED",  "efo_id": "E2",
         "phylo_score": 0.6, "n_models": 2, "sources": "impc"},
    ])


# ----------------------------------------------------------------------
# Step behaviour
# ----------------------------------------------------------------------
def test_step_max_over_drug_targets():
    """For a drug hitting multiple targets, the strongest phylo evidence wins."""
    out = steps.add_phylo_evidence(_edges(), _phylo(), _cfg())
    d1 = out[(out.drug_id == "D1") & (out.efo_id == "E1")].iloc[0]
    # D1 hits HIGH (1.0) and LOW (0.3) -> max = 1.0
    assert d1["phylo_score"] == 1.0
    assert d1["phylo_factor"] == 1.5    # 1.0 + 0.5 * 1.0
    assert d1["phylo_n_models"] == 5
    assert d1["phylo_sources"] == "impc"


def test_step_omits_drugs_without_evidence():
    """A drug with no evidence on any target is absent from the output -- merge in
    pipeline gives it NaN -> downstream defaults to factor 1.0 (no penalty)."""
    out = steps.add_phylo_evidence(_edges(), _phylo(), _cfg())
    assert not ((out.drug_id == "D2") & (out.efo_id == "E1")).any()


def test_step_boost_factor_scales():
    """factor = 1 + boost * phylo_score."""
    out = steps.add_phylo_evidence(_edges(), _phylo(), _cfg(boost=1.0))
    d3 = out[(out.drug_id == "D3") & (out.efo_id == "E2")].iloc[0]
    # MED phylo_score = 0.6, boost = 1.0 -> factor = 1.6
    assert d3["phylo_factor"] == 1.6
    out0 = steps.add_phylo_evidence(_edges(), _phylo(), _cfg(boost=0.0))
    d3 = out0[(out0.drug_id == "D3") & (out0.efo_id == "E2")].iloc[0]
    # boost=0 -> factor=1.0 even with evidence (inspection-only)
    assert d3["phylo_factor"] == 1.0


def test_step_empty_phylo_returns_empty():
    empty_phylo = pd.DataFrame(columns=["target_symbol", "efo_id", "phylo_score",
                                        "n_models", "sources"])
    out = steps.add_phylo_evidence(_edges(), empty_phylo, _cfg())
    assert len(out) == 0
    assert list(out.columns) == ["drug_id", "efo_id", "phylo_factor", "phylo_score",
                                 "phylo_n_models", "phylo_sources"]


def test_step_empty_edges_returns_empty():
    empty_edges = pd.DataFrame(columns=["drug_id", "target_symbol", "efo_id"])
    out = steps.add_phylo_evidence(empty_edges, _phylo(), _cfg())
    assert len(out) == 0


# ----------------------------------------------------------------------
# Score integration -- phylo_factor multiplies opportunity, never penalises absence
# ----------------------------------------------------------------------
def _scoring_cfg(w_inv_boost):
    """Minimal config for steps.score()."""
    return Config(
        raw={"phylogenetics": {"boost_factor": w_inv_boost},
             "direction": {}, "tissue": {}, "scoring": {}},
        root=Path("."),
    )


def _hyp(phylo_factor: float | None):
    base = {
        "drug_id": "D1", "efo_id": "E1", "disease_name": "test",
        "lead_target": "T1", "mechanistic_support": 0.8, "n_targets": 1,
        "lead_target_genecards": "", "evidence_targets": "T1(0.80)",
        "novelty": 1.0, "novelty_status": "novel",
        "direction_factor": 1.0, "direction_status": "aligned",
        "direction_alignment": 1.0,
        "tissue_factor": 1.0, "tissue_status": "expressed", "tissue_evidence": "lung",
        "n_drug_targets": 1,
    }
    if phylo_factor is not None:
        base["phylo_factor"] = phylo_factor
    return pd.DataFrame([base])


def test_score_no_phylo_column_means_no_boost():
    """When the phylo step didn't run, score() must default factor=1.0."""
    df = steps.score(_hyp(phylo_factor=None), _scoring_cfg(0.5))
    # Expected opp = 0.8 (mech) * 1.0 (novelty) * 1.0 (dir) * 1.0 (tissue) * 1.0 (phylo) / sqrt(1) = 0.8
    assert abs(df.iloc[0]["opportunity"] - 0.8) < 1e-4


def test_score_phylo_factor_boosts_opportunity():
    """phylo_factor > 1.0 multiplies opportunity."""
    df = steps.score(_hyp(phylo_factor=1.5), _scoring_cfg(0.5))
    # Expected opp = 0.8 * 1.5 = 1.2
    assert abs(df.iloc[0]["opportunity"] - 1.2) < 1e-4


def test_score_missing_phylo_factor_is_1_not_penalty():
    """NaN in phylo_factor column must fillna(1.0), NOT some unknown_factor."""
    hyp = _hyp(phylo_factor=None)
    hyp["phylo_factor"] = float("nan")
    df = steps.score(hyp, _scoring_cfg(0.5))
    assert abs(df.iloc[0]["opportunity"] - 0.8) < 1e-4


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"OK {name}")
