"""Unit tests for pathway-level mechanistic evidence (Reactome via OT)."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.config import Config  # noqa: E402
from repurpose import steps  # noqa: E402


def _cfg(enabled=True, boost=0.3, max_size=80, min_assoc=0.3, saturation=5):
    return Config(
        raw={"pathway": {"enabled": enabled, "boost_factor": boost,
                         "max_pathway_size": max_size, "min_assoc": min_assoc,
                         "saturation": saturation}},
        root=Path("."),
    )


def _ranked():
    return pd.DataFrame([
        {"rank": 1, "drug_id": "D1", "efo_id": "E1", "opportunity": 1.0},
        {"rank": 2, "drug_id": "D2", "efo_id": "E1", "opportunity": 1.0},
        {"rank": 3, "drug_id": "D3", "efo_id": "E1", "opportunity": 1.0},
    ])


def _drug_targets():
    return pd.DataFrame([
        # D1 hits T_HIT (same target as disease's strong target)
        {"drug_id": "D1", "target_symbol": "T_HIT", "is_direct": True},
        # D2 hits T_PATHWAY (different target but same pathway as T_HIT)
        {"drug_id": "D2", "target_symbol": "T_PATHWAY", "is_direct": True},
        # D3 hits T_UNRELATED (no pathway overlap)
        {"drug_id": "D3", "target_symbol": "T_UNRELATED", "is_direct": True},
    ])


def _target_pathways():
    return pd.DataFrame([
        # T_HIT and T_PATHWAY are both in P1 (small pathway, 2 members)
        {"target_symbol": "T_HIT", "pathway_id": "R-HSA-1",
         "pathway_name": "Specific pathway", "top_level": "Signaling"},
        {"target_symbol": "T_PATHWAY", "pathway_id": "R-HSA-1",
         "pathway_name": "Specific pathway", "top_level": "Signaling"},
        # T_UNRELATED is in P2 alone
        {"target_symbol": "T_UNRELATED", "pathway_id": "R-HSA-2",
         "pathway_name": "Other pathway", "top_level": "Metabolism"},
    ])


def _target_disease():
    return pd.DataFrame([
        # E1 has T_HIT as a strong target
        {"efo_id": "E1", "target_symbol": "T_HIT", "assoc_score": 0.9},
        # Weak association below min_assoc -- should NOT contribute
        {"efo_id": "E1", "target_symbol": "T_NOISE", "assoc_score": 0.1},
    ])


def test_pathway_boost_for_co_membership():
    """Drug whose target shares a pathway with disease's strong target
    gets a pathway_factor > 1.0."""
    out = steps.add_pathway_evidence(_ranked(), _drug_targets(),
                                      _target_pathways(), _target_disease(), _cfg())
    d2 = out[out.drug_id == "D2"].iloc[0]
    # D2 has 1 pathway overlap (R-HSA-1 via T_PATHWAY)
    assert d2["n_pathway_overlap"] == 1
    assert d2["pathway_factor"] > 1.0
    # Opportunity boosted
    assert d2["opportunity"] > 1.0


def test_no_boost_for_unrelated_targets():
    out = steps.add_pathway_evidence(_ranked(), _drug_targets(),
                                      _target_pathways(), _target_disease(), _cfg())
    d3 = out[out.drug_id == "D3"].iloc[0]
    assert d3["n_pathway_overlap"] == 0
    assert d3["pathway_factor"] == 1.0
    assert d3["opportunity"] == 1.0


def test_direct_hit_excluded_under_indirect_only():
    """Drug whose target IS the disease's strong target gets NO pathway boost --
    that credit already lives in mechanistic_support. Indirect-only is the
    default and the whole point of the pathway step."""
    out = steps.add_pathway_evidence(_ranked(), _drug_targets(),
                                      _target_pathways(), _target_disease(), _cfg())
    d1 = out[out.drug_id == "D1"].iloc[0]
    # D1's only pathway bridge is T_HIT <-> T_HIT (direct); excluded.
    assert d1["n_pathway_overlap"] == 0
    assert d1["pathway_factor"] == 1.0


def test_indirect_only_off_credits_direct_too():
    """With indirect_only=False, the direct-target bridge counts -- back to
    the naive co-membership behaviour."""
    cfg_naive = Config(
        raw={"pathway": {"enabled": True, "boost_factor": 0.3,
                         "max_pathway_size": 80, "min_assoc": 0.3,
                         "saturation": 5, "indirect_only": False}},
        root=Path("."),
    )
    out = steps.add_pathway_evidence(_ranked(), _drug_targets(),
                                      _target_pathways(), _target_disease(), cfg_naive)
    d1 = out[out.drug_id == "D1"].iloc[0]
    assert d1["n_pathway_overlap"] >= 1
    assert d1["pathway_factor"] > 1.0


def test_filter_excludes_generic_pathways():
    """Pathways with > max_pathway_size genes get filtered out."""
    big_pw = pd.DataFrame([
        {"target_symbol": f"GENE_{i}", "pathway_id": "R-HSA-BIG",
         "pathway_name": "Big generic", "top_level": "X"}
        for i in range(200)
    ] + [
        # Add our drug + disease genes to the big pathway too
        {"target_symbol": "T_HIT", "pathway_id": "R-HSA-BIG",
         "pathway_name": "Big generic", "top_level": "X"},
        {"target_symbol": "T_PATHWAY", "pathway_id": "R-HSA-BIG",
         "pathway_name": "Big generic", "top_level": "X"},
    ])
    out = steps.add_pathway_evidence(_ranked(), _drug_targets(), big_pw,
                                      _target_disease(),
                                      _cfg(max_size=80))
    # The big pathway is filtered; no specific pathway exists -> no overlap
    d2 = out[out.drug_id == "D2"].iloc[0]
    assert d2["n_pathway_overlap"] == 0


def test_weak_disease_associations_excluded():
    """Disease associations below min_assoc don't contribute pathways."""
    weak_td = pd.DataFrame([
        {"efo_id": "E1", "target_symbol": "T_HIT", "assoc_score": 0.05},  # weak
    ])
    out = steps.add_pathway_evidence(_ranked(), _drug_targets(),
                                      _target_pathways(), weak_td, _cfg(min_assoc=0.3))
    for _, row in out.iterrows():
        assert row["n_pathway_overlap"] == 0


def test_disabled_passes_through():
    """enabled=False -> all rows get factor 1.0, no opportunity change."""
    out = steps.add_pathway_evidence(_ranked(), _drug_targets(),
                                      _target_pathways(), _target_disease(),
                                      _cfg(enabled=False))
    assert (out["pathway_factor"] == 1.0).all()
    assert (out["opportunity"] == _ranked()["opportunity"].values).all()


def test_score_saturates_at_high_overlap():
    """Past saturation, pathway_score caps at 1.0 (factor caps at 1+boost)."""
    # 5 targets all in 5 distinct small pathways, all also have disease's strong targets
    many_pw = []
    many_td = []
    for i in range(10):
        many_pw.append({"target_symbol": "T_PATHWAY", "pathway_id": f"R-HSA-{i}",
                         "pathway_name": f"P{i}", "top_level": "Signaling"})
        many_pw.append({"target_symbol": f"DT_{i}", "pathway_id": f"R-HSA-{i}",
                         "pathway_name": f"P{i}", "top_level": "Signaling"})
        many_td.append({"efo_id": "E1", "target_symbol": f"DT_{i}", "assoc_score": 0.9})
    out = steps.add_pathway_evidence(_ranked(), _drug_targets(),
                                      pd.DataFrame(many_pw), pd.DataFrame(many_td),
                                      _cfg(saturation=3, boost=0.3))
    d2 = out[out.drug_id == "D2"].iloc[0]
    assert d2["n_pathway_overlap"] == 10
    # pathway_score should be at cap (1.0)
    assert d2["pathway_score"] == 1.0
    # pathway_factor = 1 + 0.3 * 1.0 = 1.3
    assert d2["pathway_factor"] == 1.3


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"OK {name}")
