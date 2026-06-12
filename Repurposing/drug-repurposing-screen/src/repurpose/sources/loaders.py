"""Canonical table loaders (works for both demo CSVs and full-mode files)."""
from __future__ import annotations
from pathlib import Path
import pandas as pd

_BOOL = {"true": True, "1": True, "yes": True, "y": True,
         "false": False, "0": False, "no": False, "n": False, "": False}


def _read(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _to_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().map(_BOOL).fillna(False).astype(bool)


def load_drugs(path: Path) -> pd.DataFrame:
    df = _read(path)
    df["max_phase"] = pd.to_numeric(df["max_phase"], errors="coerce").fillna(0).astype(int)
    for col in ("approved_us", "approved_eu"):
        df[col] = _to_bool(df[col])
    if "modality" not in df:
        df["modality"] = "other"
    return df[["drug_id", "drug_name", "max_phase", "approved_us", "approved_eu", "modality"]]


def load_drug_targets(path: Path) -> pd.DataFrame:
    df = _read(path)
    for col in ("action_type", "mechanism_of_action"):
        if col not in df:
            df[col] = ""
    return df[["drug_id", "target_symbol", "action_type", "mechanism_of_action"]]


def load_drug_indications(path: Path) -> pd.DataFrame:
    df = _read(path)
    if "indication_name" not in df:
        df["indication_name"] = ""
    return df[["drug_id", "efo_id", "indication_name"]]


def load_target_disease(path: Path) -> pd.DataFrame:
    df = _read(path)
    df["assoc_score"] = pd.to_numeric(df["assoc_score"], errors="coerce").fillna(0.0)
    if "disease_name" not in df:
        df["disease_name"] = ""
    return df[["target_symbol", "efo_id", "disease_name", "assoc_score"]]


def load_disease_ontology(path: Path) -> pd.DataFrame:
    df = _read(path)
    if "parent_efo_id" not in df:
        df["parent_efo_id"] = ""
    return df[["efo_id", "disease_name", "parent_efo_id"]]


def load_target_expression(path: Path) -> pd.DataFrame:
    """Baseline target expression per tissue, normalised to [0,1]. Empty if missing."""
    p = Path(path)
    if not p.exists():
        return pd.DataFrame(columns=["target_symbol", "tissue", "expression"])
    df = _read(p)
    df["expression"] = pd.to_numeric(df["expression"], errors="coerce").fillna(0.0)
    return df[["target_symbol", "tissue", "expression"]]


def load_disease_tissue(path: Path) -> pd.DataFrame:
    """Disease -> relevant tissue(s) with relevance weight in [0,1]. Empty if missing."""
    p = Path(path)
    if not p.exists():
        return pd.DataFrame(columns=["efo_id", "tissue", "relevance"])
    df = _read(p)
    df["relevance"] = pd.to_numeric(df["relevance"], errors="coerce").fillna(0.0)
    return df[["efo_id", "tissue", "relevance"]]


def load_gene_info(path: Path) -> pd.DataFrame:
    """Gene symbol -> full target name. Empty if missing."""
    p = Path(path)
    if not p.exists():
        return pd.DataFrame(columns=["symbol", "gene_name"])
    df = _read(p)
    return df[["symbol", "gene_name"]]


def load_phylo_evidence(path: Path) -> pd.DataFrame:
    """Model-organism (orthologous gene) evidence per (target, disease).

    Columns: target_symbol, efo_id, phylo_score in [0, 1], n_models,
    sources (comma-separated source datasources, e.g. "impc").
    Returns an empty frame if the file is missing -- the pipeline will
    treat all pairs as 'no phylo evidence' (factor 1.0, no penalty).
    """
    p = Path(path)
    if not p.exists():
        return pd.DataFrame(columns=["target_symbol", "efo_id", "phylo_score",
                                     "n_models", "sources"])
    df = _read(p)
    df["phylo_score"] = pd.to_numeric(df["phylo_score"], errors="coerce").fillna(0.0)
    return df[["target_symbol", "efo_id", "phylo_score", "n_models", "sources"]]


def load_target_direction(path: Path) -> pd.DataFrame:
    """Direction-of-effect: therapeutic_direction in {-1, +1}.

    +1 = increasing target activity is therapeutic (an AGONIST is wanted);
    -1 = decreasing target activity is therapeutic (an INHIBITOR is wanted).
    Returns an empty frame (no rows) if the file is missing.
    """
    p = Path(path)
    if not p.exists():
        return pd.DataFrame(columns=["target_symbol", "efo_id", "therapeutic_direction", "evidence"])
    df = _read(p)
    df["therapeutic_direction"] = pd.to_numeric(df["therapeutic_direction"], errors="coerce").fillna(0).astype(int)
    if "evidence" not in df:
        df["evidence"] = ""
    return df[["target_symbol", "efo_id", "therapeutic_direction", "evidence"]]
