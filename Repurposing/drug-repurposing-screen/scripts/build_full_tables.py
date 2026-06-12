#!/usr/bin/env python
"""Materialise canonical tables from downloaded bulk dumps (mode: full).

Edit the paths below to match your download locations, then run:
    python scripts/build_full_tables.py
Outputs data/full/*.csv in the canonical schema the pipeline expects.
"""
from __future__ import annotations
import re
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from repurpose.sources import adapters  # noqa: E402

FULL = ROOT / "data" / "full"
FULL.mkdir(parents=True, exist_ok=True)
OT_ASSOC_DIR  = FULL / "associationByOverallDirect"
OT_EVIDENCE   = FULL / "evidence"               # for direction-of-effect
OT_TARGETS    = FULL / "targets"                # for gene names + ensembl->symbol map
OT_BASELINE   = FULL / "baselineExpression"     # for tissue expression
OT_GENE_MAP   = FULL / "ot_gene_map.csv"        # ensembl_id,target_symbol
EFO_OBO       = FULL / "efo.obo"

def resolve_chembl_sqlite() -> Path:
    """Find the latest extracted ChEMBL sqlite database in data/full."""
    candidates = list(FULL.glob("chembl_*/chembl_*_sqlite/chembl_*.db"))
    if not candidates:
        raise FileNotFoundError(
            "No ChEMBL sqlite DB found under data/full/chembl_*/chembl_*_sqlite/chembl_*.db"
        )
    def version_key(path: Path) -> int:
        m = re.search(r"chembl_(\d+)\.db$", path.name)
        return int(m.group(1)) if m else -1
    return max(candidates, key=version_key)


def ensure_ot_gene_map(targets_dir: Path, out_csv: Path) -> Path:
    """Create ensembl_id -> target_symbol map from OT targets parquet if missing."""
    if out_csv.exists():
        return out_csv
    targets = pd.read_parquet(parquet_paths(targets_dir), columns=["id", "approvedSymbol"])
    gene_map = (targets.rename(columns={"id": "ensembl_id", "approvedSymbol": "target_symbol"})
                      [["ensembl_id", "target_symbol"]]
                      .dropna()
                      .drop_duplicates())
    gene_map.to_csv(out_csv, index=False)
    return out_csv


def parquet_paths(dataset_dir: Path) -> list[str]:
    """Return recursive parquet shard paths for a downloaded OT dataset."""
    paths = sorted(dataset_dir.rglob("*.parquet"))
    if not paths:
        raise FileNotFoundError(f"No parquet files found under {dataset_dir}")
    return [str(p) for p in paths]

def export_filtered_csv(df: pd.DataFrame, out_path: Path, key_cols: list[str], label: str) -> None:
    """Drop rows with null values in key columns, then export to CSV."""
    missing = [c for c in key_cols if c not in df.columns]
    if missing:
        raise KeyError(f"{label}: missing required key columns: {missing}")
    filtered = df.copy()
    null_like_strings = {"", "NA", "N/A", "NULL", "NONE", "NAN"}
    for col in key_cols:
        if pd.api.types.is_object_dtype(filtered[col]) or pd.api.types.is_string_dtype(filtered[col]):
            normalized = filtered[col].astype("string").str.strip()
            normalized = normalized.mask(normalized.str.upper().isin(null_like_strings), pd.NA)
            filtered[col] = normalized
    dropped = int(filtered[key_cols].isna().any(axis=1).sum())
    if dropped:
        print(f"{label}: dropping {dropped} rows with null key values in {key_cols}")
    filtered.dropna(subset=key_cols).to_csv(out_path, index=False)


def main() -> None:
    chembl_sqlite = resolve_chembl_sqlite()
    gene_map_csv = ensure_ot_gene_map(OT_TARGETS, OT_GENE_MAP)
    ot_assoc = parquet_paths(OT_ASSOC_DIR)
    ot_targets = parquet_paths(OT_TARGETS)
    ot_evidence = parquet_paths(OT_EVIDENCE)
    ot_baseline = parquet_paths(OT_BASELINE)
    print(f"Using ChEMBL DB: {chembl_sqlite}")
    print(f"Using OT gene map: {gene_map_csv}")
    print("ChEMBL -> drugs.csv")
    export_filtered_csv(
        adapters.chembl_drugs(str(chembl_sqlite)),
        FULL / "drugs.csv",
        ["drug_id", "drug_name"],
        "drugs.csv",
    )
    print("ChEMBL -> drug_targets.csv")
    export_filtered_csv(
        adapters.chembl_drug_targets(str(chembl_sqlite)),
        FULL / "drug_targets.csv",
        ["drug_id", "target_symbol"],
        "drug_targets.csv",
    )
    print("ChEMBL -> drug_indications.csv")
    export_filtered_csv(
        adapters.chembl_drug_indications(str(chembl_sqlite)),
        FULL / "drug_indications.csv",
        ["drug_id", "efo_id"],
        "drug_indications.csv",
    )
    print("Open Targets -> target_disease.csv")
    export_filtered_csv(
        adapters.opentargets_target_disease(ot_assoc, str(gene_map_csv)),
        FULL / "target_disease.csv",
        ["target_symbol", "efo_id", "assoc_score"],
        "target_disease.csv",
    )
    print("Open Targets -> target_direction.csv")
    export_filtered_csv(
        adapters.opentargets_target_direction(ot_evidence, str(gene_map_csv)),
        FULL / "target_direction.csv",
        ["target_symbol", "efo_id", "therapeutic_direction"],
        "target_direction.csv",
    )
    print("Open Targets -> gene_info.csv")
    export_filtered_csv(
        adapters.opentargets_gene_info(ot_targets),
        FULL / "gene_info.csv",
        ["symbol", "gene_name"],
        "gene_info.csv",
    )
    print("Open Targets -> target_expression.csv")
    export_filtered_csv(
        adapters.opentargets_target_expression(ot_baseline, str(gene_map_csv)),
        FULL / "target_expression.csv",
        ["target_symbol", "tissue", "expression"],
        "target_expression.csv",
    )
    print("EFO -> disease_ontology.csv")
    export_filtered_csv(
        adapters.efo_ontology(str(EFO_OBO)),
        FULL / "disease_ontology.csv",
        ["efo_id", "parent_efo_id"],
        "disease_ontology.csv",
    )
    print("NOTE: disease_tissue.csv (disease -> UBERON tissue) is curated/derived; see README.")
    print("Done. Set mode: full in config.yaml and point paths at data/full/*.")


if __name__ == "__main__":
    main()
