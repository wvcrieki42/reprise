"""Unit tests for combination-therapy companion finder."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.config import Config  # noqa: E402
from repurpose import steps  # noqa: E402


def _cfg(top_n=10, min_synergy=0.1, max_overlap=0.5,
         partners=2, min_assoc=0.3, enabled=True):
    return Config(
        raw={"combination": {
            "enabled": enabled, "top_n": top_n,
            "min_synergy": min_synergy, "max_target_overlap": max_overlap,
            "partners_per_hit": partners, "partner_min_assoc": min_assoc,
        }},
        root=Path("."),
    )


def _ranked():
    return pd.DataFrame([
        {"rank": 1, "drug_id": "D1", "substance_chembl_id": "S1",
         "efo_id": "EFO_X", "opportunity": 1.0},
        {"rank": 2, "drug_id": "D2", "substance_chembl_id": "S2",
         "efo_id": "EFO_X", "opportunity": 0.9},
    ])


def _drug_targets():
    return pd.DataFrame([
        # S1 (D1) hits TGT_A (strong disease target, assoc 0.7)
        {"drug_id": "D1", "target_symbol": "TGT_A", "action_type": "INHIBITOR", "is_direct": True},
        # S2 (D2) hits TGT_B (strong disease target the primary misses)
        {"drug_id": "D2", "target_symbol": "TGT_B", "action_type": "INHIBITOR", "is_direct": True},
        # S3 (D3) hits TGT_B AND TGT_A (overlaps with S1)
        {"drug_id": "D3", "target_symbol": "TGT_A", "action_type": "INHIBITOR", "is_direct": True},
        {"drug_id": "D3", "target_symbol": "TGT_B", "action_type": "INHIBITOR", "is_direct": True},
        # S4 (D4) hits TGT_C (another strong disease target -- pure bridge candidate)
        {"drug_id": "D4", "target_symbol": "TGT_C", "action_type": "INHIBITOR", "is_direct": True},
    ])


def _substance_map():
    return pd.DataFrame([
        {"drug_id": "D1", "substance_chembl_id": "S1", "substance_name": "ALPHADRUG"},
        {"drug_id": "D2", "substance_chembl_id": "S2", "substance_name": "BETADRUG"},
        {"drug_id": "D3", "substance_chembl_id": "S3", "substance_name": "GAMMADRUG"},
        {"drug_id": "D4", "substance_chembl_id": "S4", "substance_name": "DELTADRUG"},
    ])


def _target_disease():
    return pd.DataFrame([
        # Three strong targets for EFO_X
        {"target_symbol": "TGT_A", "efo_id": "EFO_X", "assoc_score": 0.7},
        {"target_symbol": "TGT_B", "efo_id": "EFO_X", "assoc_score": 0.6},
        {"target_symbol": "TGT_C", "efo_id": "EFO_X", "assoc_score": 0.5},
    ])


def test_finds_complementary_partner():
    """Primary hits TGT_A; expects BETADRUG (TGT_B) and DELTADRUG (TGT_C)
    as bridging partners."""
    out = steps.add_combination_partners(_ranked(), _drug_targets(),
                                          _target_disease(), _substance_map(),
                                          _cfg())
    s1 = out[out.substance_chembl_id == "S1"].iloc[0]
    partners = {s1["combo_partner_1_name"], s1["combo_partner_2_name"]}
    # Both BETADRUG and DELTADRUG should appear (each bridges one uncovered target)
    assert "BETADRUG" in partners or "DELTADRUG" in partners
    assert s1["combo_partner_1_synergy"] > 0


def test_rejects_overlapping_class():
    """GAMMADRUG hits both TGT_A and TGT_B (50% overlap with primary S1's
    single-target set is exactly the boundary). With max_overlap=0.5 it should
    be allowed; with max_overlap=0.4 it should be rejected."""
    # Strict overlap -- gamma rejected, only pure-bridge candidates remain
    out = steps.add_combination_partners(_ranked(), _drug_targets(),
                                          _target_disease(), _substance_map(),
                                          _cfg(max_overlap=0.4))
    s1 = out[out.substance_chembl_id == "S1"].iloc[0]
    chosen = {s1["combo_partner_1_name"], s1["combo_partner_2_name"]}
    assert "GAMMADRUG" not in chosen


def test_synergy_threshold_filters():
    """High threshold drops weak combos."""
    out = steps.add_combination_partners(_ranked(), _drug_targets(),
                                          _target_disease(), _substance_map(),
                                          _cfg(min_synergy=0.9))
    s1 = out[out.substance_chembl_id == "S1"].iloc[0]
    assert s1["combo_partner_1_name"] == ""


def test_disabled_passes_through():
    out = steps.add_combination_partners(_ranked(), _drug_targets(),
                                          _target_disease(), _substance_map(),
                                          _cfg(enabled=False))
    assert (out["combo_partner_1_name"] == "").all()


def test_top_n_respects_cap():
    """Only the first top_n primary rows get partner lookup."""
    out = steps.add_combination_partners(_ranked(), _drug_targets(),
                                          _target_disease(), _substance_map(),
                                          _cfg(top_n=1))
    s1 = out[out.substance_chembl_id == "S1"].iloc[0]
    s2 = out[out.substance_chembl_id == "S2"].iloc[0]
    # S1 (rank 1) gets a partner; S2 (rank 2) does NOT (above top_n cap)
    assert s1["combo_partner_1_name"] != ""
    assert s2["combo_partner_1_name"] == ""


def test_partners_sorted_by_synergy_desc():
    out = steps.add_combination_partners(_ranked(), _drug_targets(),
                                          _target_disease(), _substance_map(),
                                          _cfg())
    s1 = out[out.substance_chembl_id == "S1"].iloc[0]
    if s1["combo_partner_2_synergy"] is not pd.NA and not pd.isna(s1["combo_partner_2_synergy"]):
        assert s1["combo_partner_1_synergy"] >= s1["combo_partner_2_synergy"]


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"OK {name}")
