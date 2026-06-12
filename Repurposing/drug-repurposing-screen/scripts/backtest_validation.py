"""Backtest the screen against known drug-repurposing successes.

For each curated case (data/curated/repurposing_validation.yaml), look up
whether the pipeline's mechanism layer would have flagged the (drug,
repurposed_indication) connection independent of novelty subtraction.

Specifically, for each case we:
  1. Find the drug in drugs.csv (ChEMBL pref_name exact match, fallback to
     case-insensitive contains).
  2. Get the drug's direct + STRING-expanded targets from drug_targets.csv.
  3. Resolve the repurposed disease to an OT efo_id -- explicit `efo`
     wins; otherwise a case-insensitive name match against
     target_disease.csv.
  4. For each (target, disease) edge, pull assoc_score from
     target_disease.csv.
  5. Compute mech_support exactly as the pipeline does: noisy-OR of
     `target_weight * assoc_score` over the drug's targets.

The metric we care about is "would the screen surface this?". We define
HIT as mech_support >= 0.3 (matches the typical opportunity threshold
once novelty etc. are 1.0). Lower hits explain why the screen would
miss them; higher hits are wins.

We also report the per-target contribution breakdown so failures are
debuggable -- "drug not found", "disease not found", "no informative
target overlap", etc.

Usage:
    python scripts/backtest_validation.py
        [--config config.full.yaml]
        [--cases data/curated/repurposing_validation.yaml]
        [--threshold 0.3]
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from repurpose.config import load_config  # noqa: E402
from repurpose.sources import loaders  # noqa: E402


def _find_drug_ids(name: str, drugs: pd.DataFrame) -> list[str]:
    """Return ALL ChEMBL IDs whose drug_name matches the input -- parent + salts.

    ChEMBL's mechanism_of_action coverage is patchy at the parent level
    (SILDENAFIL CHEMBL192 has no MoA rows; SILDENAFIL CITRATE CHEMBL1737
    does), so we union targets across every formulation that shares a
    canonical name. Matches by uppercase substring against drug_name.
    """
    name_up = name.strip().upper()
    exact = drugs[drugs["drug_name"].str.upper() == name_up]
    # Always also include salt-form variants ("SILDENAFIL CITRATE",
    # "MEMANTINE HYDROCHLORIDE" etc.) -- the MoA / target rows often live
    # on the salt, not the parent.
    contains = drugs[drugs["drug_name"].str.upper().str.contains(
        name_up, na=False, regex=False)]
    return (pd.concat([exact, contains])
              .drop_duplicates("drug_id")["drug_id"].tolist())


def _resolve_efo(case: dict, target_disease: pd.DataFrame) -> tuple[str, str]:
    """Returns (efo_id, matched_disease_name) or ('', '')."""
    if case.get("repurposed_disease_efo"):
        efo = case["repurposed_disease_efo"]
        rows = target_disease[target_disease["efo_id"] == efo]
        if len(rows):
            return efo, rows.iloc[0]["disease_name"]
    name = (case.get("repurposed_disease_name") or "").strip().lower()
    if not name:
        return "", ""
    # Exact-only -- the previous substring fallback yielded false matches
    # like "Antimigraine preparation use measurement" for migraine.
    td = target_disease[target_disease["disease_name"].str.lower() == name]
    if len(td):
        return td.iloc[0]["efo_id"], td.iloc[0]["disease_name"]
    return "", ""


def _mech_support_for(drug_ids: list[str], efo_id: str,
                       drug_targets: pd.DataFrame,
                       target_disease: pd.DataFrame,
                       min_assoc: float = 0.1) -> dict:
    """Mirror steps.disease_edges + propagate_disease(noisy_or) for one pair."""
    dt = drug_targets[drug_targets["drug_id"].isin(drug_ids)].copy()
    if dt.empty:
        return {"mech_support": 0.0, "n_targets": 0, "edges": []}
    td = target_disease[(target_disease["efo_id"] == efo_id)
                        & (target_disease["assoc_score"] >= min_assoc)].copy()
    if td.empty:
        return {"mech_support": 0.0, "n_targets": 0, "edges": []}
    # The pipeline weights direct targets at 1.0 and STRING neighbours at
    # neighbor_weight * string_score; for backtest validity we use only
    # direct targets (string expansion is a separate, optional layer).
    if "is_direct" in dt.columns:
        dt = dt[dt["is_direct"].fillna(True)]
        if dt.empty:
            return {"mech_support": 0.0, "n_targets": 0, "edges": []}
    if "target_weight" not in dt.columns:
        dt["target_weight"] = 1.0
    e = dt.merge(td, on="target_symbol", how="inner")
    if e.empty:
        return {"mech_support": 0.0, "n_targets": 0, "edges": []}
    e["contrib"] = (e["target_weight"].clip(0, 1).astype(float)
                    * e["assoc_score"].clip(0, 1).astype(float))
    mech_support = float(1.0 - np.prod(1.0 - np.clip(e["contrib"], 0, 0.999)))
    edges = (e.sort_values("contrib", ascending=False)
              [["target_symbol", "assoc_score", "contrib"]]
              .head(8).to_dict("records"))
    return {"mech_support": round(mech_support, 4),
            "n_targets": int(e["target_symbol"].nunique()),
            "edges": edges}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.full.yaml")
    ap.add_argument("--cases",
                    default="data/curated/repurposing_validation.yaml")
    ap.add_argument("--threshold", type=float, default=0.3)
    args = ap.parse_args()

    cfg = load_config(ROOT / args.config)
    drugs = loaders.load_drugs(cfg.path("drugs"))
    drug_targets = loaders.load_drug_targets(cfg.path("drug_targets"))
    target_disease = loaders.load_target_disease(cfg.path("target_disease"))
    print(f"loaded: {len(drugs):,} drugs, {len(drug_targets):,} drug-target "
          f"edges, {len(target_disease):,} target-disease associations")
    print()

    cases = yaml.safe_load((ROOT / args.cases).read_text())["cases"]
    print(f"validating {len(cases)} curated repurposing cases at threshold "
          f"mech_support >= {args.threshold}")
    print()

    results = []
    for case in cases:
        drug_ids = _find_drug_ids(case["drug"], drugs)
        if not drug_ids:
            results.append({"drug": case["drug"],
                            "disease": case.get("repurposed_disease_name", ""),
                            "status": "DRUG_NOT_IN_UNIVERSE",
                            "mech_support": 0.0, "n_targets": 0, "edges": []})
            continue
        efo_id, matched_name = _resolve_efo(case, target_disease)
        if not efo_id:
            results.append({"drug": case["drug"], "drug_id": ",".join(drug_ids),
                            "disease": case.get("repurposed_disease_name", ""),
                            "status": "DISEASE_NOT_IN_OT",
                            "mech_support": 0.0, "n_targets": 0, "edges": []})
            continue
        mech = _mech_support_for(drug_ids, efo_id,
                                   drug_targets, target_disease)
        expected = set(case.get("mechanism_targets") or [])
        actual = {e["target_symbol"] for e in mech["edges"]}
        hit = mech["mech_support"] >= args.threshold
        results.append({
            "drug": case["drug"], "drug_id": ",".join(drug_ids),
            "disease": matched_name, "efo": efo_id,
            "expected_targets": ",".join(sorted(expected)),
            "actual_top_targets": ",".join(sorted(actual)),
            "expected_overlap": bool(expected & actual),
            "mech_support": mech["mech_support"],
            "n_targets": mech["n_targets"],
            "status": "HIT" if hit else "MISS",
            "edges": mech["edges"],
            "year": case.get("year_repurposed"),
        })

    # Per-case report
    print(f"{'STATUS':<22} {'DRUG':<28} {'->':>4} {'DISEASE':<40} {'mech':>7} {'overlap':>8}")
    print("-" * 120)
    for r in results:
        print(f"{r['status']:<22} {r['drug'][:27]:<28} {'->':>4} "
              f"{(r.get('disease') or '')[:39]:<40} "
              f"{r.get('mech_support', 0):>7.3f} "
              f"{('YES' if r.get('expected_overlap') else '-'):>8}")
    print()

    # Summary
    statuses = pd.Series([r["status"] for r in results]).value_counts()
    print("=== Summary ===")
    for st, n in statuses.items():
        print(f"  {st:<22} {n:>3}")
    pass_rate = statuses.get("HIT", 0) / len(results) * 100
    print(f"\n  Hit rate: {statuses.get('HIT', 0)}/{len(results)}  "
          f"({pass_rate:.0f}%)")
    print()

    # Detailed misses for debugging
    misses = [r for r in results if r["status"] == "MISS"]
    if misses:
        print("=== Cases that scored below threshold ===")
        for r in misses:
            print(f"  {r['drug']} -> {r['disease']}")
            print(f"     mech_support={r['mech_support']:.3f}, "
                  f"n_targets={r['n_targets']}")
            for e in r["edges"][:3]:
                print(f"      {e['target_symbol']:>8s}  assoc={e['assoc_score']:.3f}  contrib={e['contrib']:.3f}")
        print()

    # Cases where the top mechanism targets DON'T overlap with the curated
    # expectation -- these are interesting because the screen surfaces the
    # connection via a different bridge than the literature attributes.
    surprising = [r for r in results
                  if r["status"] == "HIT" and not r.get("expected_overlap")]
    if surprising:
        print("=== Surprising hits (high mech via different targets than the "
              "literature attributes) ===")
        for r in surprising:
            print(f"  {r['drug']} -> {r['disease']}")
            print(f"     expected: {r['expected_targets']}")
            print(f"     actually surfaced via: {r['actual_top_targets']}")
        print()


if __name__ == "__main__":
    main()
