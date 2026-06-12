"""End-to-end orchestration of the repurposing screen."""
from __future__ import annotations
import json
import time
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

from .config import Config, load_config
from .sources import loaders
from . import steps

GENECARDS = "https://www.genecards.org/cgi-bin/carddisp.pl?gene="

OUTPUT_COLS = ["rank", "drug_id", "drug_name", "substance_chembl_id", "substance_name",
               "n_variants", "variant_names",
               "modality", "efo_id", "disease_name",
               "lead_target", "lead_target_name", "lead_target_genecards",
               "mechanistic_support", "novelty", "novelty_status",
               "direction_factor", "direction_status",
               "tissue_factor", "tissue_status", "tissue_evidence",
               "phylo_factor", "phylo_score", "phylo_n_models", "phylo_sources",
               "disease_gene_match",
               "pubmed_count", "europepmc_count", "trial_count", "patent_count",
               "investigation_prior",
               "us_patients", "us_prevalence_per_100k", "is_orphan",
               "market_source", "as_of",
               "n_targets", "n_drug_targets", "opportunity",
               "evidence_targets", "data_version"]


def _maybe_literature_pass(ranked: pd.DataFrame, cfg: Config, prep: dict, log) -> pd.DataFrame:
    """Construct a LiteraturePriorClient from config and run the two-pass step.

    Enriches the top-N ranked hypotheses with a `target_synonyms` column built
    from `gene_info.gene_name` so the PubMed / Europe PMC queries OR the full
    gene name in alongside the symbol -- materially better recall for genes
    with descriptive names ("EGFR" OR "epidermal growth factor receptor").
    """
    if not cfg.get("literature", "enabled", default=False):
        return ranked
    from .sources.literature import LiteraturePriorClient
    cache_dir = cfg.root / cfg.get("literature", "cache_dir", default=".cache/literature")
    lit = LiteraturePriorClient(
        cache_dir=cache_dir,
        ncbi_api_key=cfg.get("literature", "ncbi_api_key", default=None),
        lens_api_token=cfg.get("literature", "lens_api_token", default=None),
        enable_pubmed=cfg.get("literature", "enable_pubmed", default=True),
        enable_europepmc=cfg.get("literature", "enable_europepmc", default=True),
        enable_clinicaltrials=cfg.get("literature", "enable_clinicaltrials", default=True),
        enable_lens=cfg.get("literature", "enable_lens", default=False),
        max_synonyms_target=cfg.get("literature", "max_synonyms_target", default=4),
        max_synonyms_disease=cfg.get("literature", "max_synonyms_disease", default=4),
        max_workers=cfg.get("literature", "max_workers", default=4),
    )
    n_query = min(int(cfg.get("literature", "top_n", default=5000)), len(ranked))
    log(f"literature pass: querying top {n_query} of {len(ranked)} hypotheses "
        f"(cache_dir={cache_dir}; {lit.credentials_summary()})")
    # Attach gene_name as a target synonym so the client's OR expansion picks it up.
    gi = prep.get("gene_info")
    if gi is not None and {"symbol", "gene_name"}.issubset(gi.columns):
        syn = gi[["symbol", "gene_name"]].dropna()
        syn = syn[syn["gene_name"].str.strip() != ""].rename(
            columns={"symbol": "lead_target", "gene_name": "target_synonyms"})
        ranked = ranked.merge(syn, on="lead_target", how="left")
    # Same on the disease side: OT exact synonyms keyed by efo_id.
    dsyn = prep.get("disease_synonyms")
    if dsyn is not None and not dsyn.empty:
        ranked = ranked.merge(dsyn, on="efo_id", how="left")
        ranked["disease_synonyms"] = ranked["disease_synonyms"].fillna("")
    return steps.add_literature_pass(ranked, lit, cfg)


def _maybe_market_pass(ranked: pd.DataFrame, cfg: Config, log) -> pd.DataFrame:
    """Attach US market data (patient count, prevalence, orphan flag) from a curated CSV.

    Inspection-only: this pass never modifies opportunity or rank -- it
    just adds columns the caller can sort/filter by downstream.
    """
    if not cfg.get("market", "enabled", default=False):
        return ranked
    from .sources.market import MarketSizeClient
    curated = cfg.path("disease_prevalence") if cfg.get("paths", "disease_prevalence", default=None) else None
    if curated is None:
        curated = cfg.root / cfg.get("market", "curated_csv", default="data/curated/disease_prevalence.csv")
    client = MarketSizeClient(curated_csv=curated)
    n_query = min(int(cfg.get("market", "top_n", default=10_000)), len(ranked))
    log(f"market pass: looking up US market size for top {n_query} of {len(ranked)} "
        f"hypotheses (csv={curated})")
    return steps.add_market_size(ranked, client, cfg)


def _logger(verbose: bool):
    return (lambda m: print(f"[{time.strftime('%H:%M:%S')}] {m}")) if verbose else (lambda m: None)


def _data_version(cfg: Config) -> str:
    p = cfg.get("provenance", default={}) or {}
    return (f"ChEMBL {p.get('chembl_version', '?')}; "
            f"OpenTargets {p.get('opentargets_version', '?')}; "
            f"STRING {p.get('string_version', '?')}; "
            f"EFO {p.get('efo_version', '?')}")


def _prepare(cfg: Config, log):
    """Shared preparation used by both engines (everything except the big join)."""
    drugs = loaders.load_drugs(cfg.path("drugs"))
    drug_targets_raw = loaders.load_drug_targets(cfg.path("drug_targets"))
    indications = loaders.load_drug_indications(cfg.path("drug_indications"))
    ontology = loaders.load_disease_ontology(cfg.path("disease_ontology"))

    direction_on = bool(cfg.get("direction", "enabled", default=False))
    target_direction = (loaders.load_target_direction(cfg.path("target_direction"))
                        if direction_on else
                        pd.DataFrame(columns=["target_symbol", "efo_id", "therapeutic_direction"]))
    tissue_on = bool(cfg.get("tissue", "enabled", default=False))
    target_expression = loaders.load_target_expression(cfg.path("target_expression")) if tissue_on \
        else pd.DataFrame(columns=["target_symbol", "tissue", "expression"])
    disease_tissue = loaders.load_disease_tissue(cfg.path("disease_tissue")) if tissue_on \
        else pd.DataFrame(columns=["efo_id", "tissue", "relevance"])
    phylo_on = bool(cfg.get("phylogenetics", "enabled", default=False))
    phylo_evidence = (loaders.load_phylo_evidence(cfg.path("phylo_evidence")) if phylo_on
                      else pd.DataFrame(columns=["target_symbol", "efo_id", "phylo_score",
                                                 "n_models", "sources"]))
    gene_info = loaders.load_gene_info(cfg.path("gene_info"))
    substance_map = (loaders.load_substance_map(cfg.path("chembl_substance_map"))
                     if cfg.get("paths", "chembl_substance_map", default=None) is not None
                     else pd.DataFrame(columns=["drug_id", "substance_chembl_id", "substance_name"]))
    disease_synonyms = (loaders.load_disease_synonyms(cfg.path("disease_synonyms"))
                        if cfg.get("paths", "disease_synonyms", default=None) is not None
                        else pd.DataFrame(columns=["efo_id", "disease_synonyms"]))

    universe = steps.build_universe(drugs, cfg)
    log(f"universe (approved US/EU): {len(universe)} drugs")
    dt = steps.build_drug_targets(universe, drug_targets_raw, cfg)
    log(f"drug-target edges after expansion: {len(dt)} "
        f"({int(dt['is_direct'].sum())} direct, {int((~dt['is_direct']).sum())} neighbour)")

    # Prune target_expression to only (tissues, targets) the join can actually reach.
    # At full scale (53k disease_tissue rows + 4.8M expression rows) the unfiltered
    # cross-join blows DuckDB's temp budget; restricting to the ~32 tissues in
    # disease_tissue and the ~15k targets in the drug-target edges cuts it by >10x.
    if tissue_on and not target_expression.empty:
        rel_tissues = set(disease_tissue["tissue"].dropna())
        rel_targets = set(dt["target_symbol"].dropna())
        target_expression = target_expression[
            target_expression["tissue"].isin(rel_tissues)
            & target_expression["target_symbol"].isin(rel_targets)
        ].reset_index(drop=True)
        log(f"target_expression pruned to relevant (tissue, target) pairs: "
            f"{len(target_expression):,} rows")

    indications = steps.expand_indications_by_target_class(
        indications, drug_targets_raw, cfg)
    known_exp = steps.known_expanded(indications, ontology, cfg)
    direct = drug_targets_raw.merge(universe[["drug_id"]], on="drug_id", how="inner")
    breadth = (direct.groupby("drug_id")["target_symbol"].nunique()
                     .rename("n_drug_targets").reset_index())
    return {"universe": universe, "drug_targets": dt, "target_direction": target_direction,
            "target_expression": target_expression, "disease_tissue": disease_tissue,
            "phylo_evidence": phylo_evidence, "substance_map": substance_map,
            "disease_synonyms": disease_synonyms,
            "gene_info": gene_info, "known_exp": known_exp, "breadth": breadth,
            "direction_on": direction_on, "tissue_on": tissue_on, "phylo_on": phylo_on}


def run(cfg: Config, *, verbose: bool = True) -> pd.DataFrame:
    log = _logger(verbose)
    engine = cfg.get("engine", default="pandas")
    prep = _prepare(cfg, log)

    if engine == "duckdb":
        from .backends.duckdb_engine import run_duckdb
        log("engine: duckdb (streaming target-disease from disk)")
        ranked = run_duckdb(cfg, prep, log)
    else:
        log("engine: pandas (in-memory)")
        target_disease = loaders.load_target_disease(cfg.path("target_disease"))
        edges = steps.disease_edges(prep["drug_targets"], target_disease, cfg)
        log(f"drug x target x disease edges: {len(edges)}")
        hyp = steps.propagate_disease(edges, cfg)
        hyp = steps.add_novelty(hyp, prep["known_exp"])
        if prep["direction_on"]:
            hyp = hyp.merge(steps.add_direction(edges, prep["target_direction"], cfg),
                            on=["drug_id", "efo_id"], how="left")
        if prep["tissue_on"]:
            hyp = hyp.merge(steps.add_tissue(edges, prep["target_expression"],
                                             prep["disease_tissue"], cfg),
                            on=["drug_id", "efo_id"], how="left")
        if prep["phylo_on"]:
            hyp = hyp.merge(steps.add_phylo_evidence(edges, prep["phylo_evidence"], cfg),
                            on=["drug_id", "efo_id"], how="left")
        hyp = hyp.merge(prep["breadth"], on="drug_id", how="left")
        hyp["n_drug_targets"] = hyp["n_drug_targets"].fillna(1).astype(int)
        ranked = steps.score(hyp, cfg)

    # Bring drug_name forward so the substance collapse can report variant
    # names instead of bare ChEMBL ids in `variant_names`. (The duckdb engine
    # builds ranked without drug_name; the merge in _finalize otherwise comes
    # too late.)
    if "drug_name" not in ranked.columns:
        ranked = ranked.merge(prep["universe"][["drug_id", "drug_name"]],
                              on="drug_id", how="left")

    # Collapse formulation variants up to their active-ingredient parent BEFORE
    # the downstream passes so the literature / market lookups operate on
    # substances (avoids querying the same (lead_target, disease) for 11
    # insulin variants).
    pre_collapse = len(ranked)
    ranked = steps.collapse_to_substances(ranked, prep["substance_map"], cfg)
    if len(ranked) != pre_collapse:
        log(f"substance collapse: {pre_collapse} -> {len(ranked)} rows "
            f"(grouped formulation variants)")
    # Flag (drug, disease) pairs where one of the drug's direct targets is named
    # in the disease itself -- 'drug targets the broken receptor' class. Always
    # on, inspection-only (no automatic score change).
    direct_dt = prep["drug_targets"][prep["drug_targets"].get("is_direct", True) == True] \
        if "is_direct" in prep["drug_targets"].columns else prep["drug_targets"]
    concord = steps.flag_disease_gene_concordance(
        ranked, direct_dt, prep["gene_info"])
    ranked = ranked.merge(concord, on=["drug_id", "efo_id"], how="left")
    ranked["disease_gene_match"] = ranked["disease_gene_match"].fillna("")
    ranked = _maybe_literature_pass(ranked, cfg, prep, log)
    ranked = _maybe_market_pass(ranked, cfg, log)
    ranked = _finalize(ranked, cfg, prep)
    log(f"final ranked hypotheses: {len(ranked)}")
    return ranked


def _finalize(ranked: pd.DataFrame, cfg: Config, prep: dict) -> pd.DataFrame:
    # Don't clobber a drug_name that has already been replaced by substance_name
    # during the collapse pass; only pull in modality (and drug_name if it's
    # still missing).
    cols = ["drug_id"] + [c for c in ["drug_name", "modality"]
                          if c not in ranked.columns]
    if len(cols) > 1:
        ranked = ranked.merge(prep["universe"][cols], on="drug_id", how="left")
    else:
        ranked = ranked.merge(prep["universe"][["drug_id", "modality"]],
                              on="drug_id", how="left")
    # lead-target name + GeneCards link
    gi = prep["gene_info"].rename(columns={"symbol": "lead_target", "gene_name": "lead_target_name"})
    ranked = ranked.merge(gi, on="lead_target", how="left")
    ranked["lead_target_name"] = ranked["lead_target_name"].fillna("")
    ranked["lead_target_genecards"] = ranked["lead_target"].apply(
        lambda s: (GENECARDS + str(s)) if pd.notna(s) and str(s) else "")
    ranked["data_version"] = _data_version(cfg)
    for c in OUTPUT_COLS:
        if c not in ranked.columns:
            ranked[c] = pd.NA
    return ranked[OUTPUT_COLS]


def run_from_file(config_path: str | Path, verbose: bool = True) -> pd.DataFrame:
    cfg = load_config(config_path)
    ranked = run(cfg, verbose=verbose)
    out_csv = cfg.outdir / "repurposing_hypotheses.csv"
    ranked.to_csv(out_csv, index=False)
    try:
        ranked.to_parquet(cfg.outdir / "repurposing_hypotheses.parquet", index=False)
    except Exception:
        pass
    meta = {
        "run_date_utc": datetime.now(timezone.utc).isoformat(),
        "mode": cfg.mode,
        "engine": cfg.get("engine", default="pandas"),
        "data_version": _data_version(cfg),
        "provenance": cfg.get("provenance", default={}),
        "n_hypotheses": int(len(ranked)),
        "thresholds": {
            "min_assoc": cfg.get("propagation", "min_assoc", default=0.1),
            "aggregate": cfg.get("propagation", "aggregate", default="noisy_or"),
            "novelty_radius": cfg.get("novelty", "ontology_radius", default=1),
            "direction_enabled": cfg.get("direction", "enabled", default=False),
            "tissue_enabled": cfg.get("tissue", "enabled", default=False),
            "phylogenetics_enabled": cfg.get("phylogenetics", "enabled", default=False),
            "literature_enabled": cfg.get("literature", "enabled", default=False),
            "literature_top_n": cfg.get("literature", "top_n", default=5000),
            "w_investigation": cfg.get("scoring", "w_investigation", default=0.0),
            "market_enabled": cfg.get("market", "enabled", default=False),
            "rare_disease_us_threshold": cfg.get(
                "market", "rare_disease_us_threshold", default=200_000),
            "min_opportunity": cfg.get("scoring", "min_opportunity", default=0.0),
        },
    }
    (cfg.outdir / "run_metadata.json").write_text(json.dumps(meta, indent=2, default=str))
    if verbose:
        print(f"[done] wrote {out_csv} and run_metadata.json")
    return ranked
