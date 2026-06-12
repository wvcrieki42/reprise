"""Pipeline steps. Each is a pure function over pandas DataFrames.

The expensive object is `disease_edges` (drug x target x disease). Both the
mechanistic-support aggregation and the directionality model are computed from
it, so the two engines (pandas / duckdb) stay consistent.
"""
from __future__ import annotations
import re
import numpy as np
import pandas as pd

from .config import Config
from .ontology import DiseaseOntology
from .sources.string_api import StringClient


# ----------------------------------------------------------------------
# Step 1 - approved-drug universe
# ----------------------------------------------------------------------
def build_universe(drugs: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    req_phase = cfg.get("universe", "require_max_phase", default=4)
    regions = cfg.get("universe", "regions", default=["us", "eu"])
    modalities = set(cfg.get("universe", "modalities", default=[]))
    mask = drugs["max_phase"] >= req_phase
    region_mask = np.zeros(len(drugs), dtype=bool)
    if "us" in regions:
        region_mask |= drugs["approved_us"].values
    if "eu" in regions:
        region_mask |= drugs["approved_eu"].values
    mask &= region_mask
    if modalities:
        mask &= drugs["modality"].isin(modalities).values
    return drugs[mask].copy().reset_index(drop=True)


# ----------------------------------------------------------------------
# Step 2 - drug -> targets, optionally expanded by STRING network
# ----------------------------------------------------------------------
def build_drug_targets(universe: pd.DataFrame, drug_targets: pd.DataFrame,
                       cfg: Config) -> pd.DataFrame:
    """Return drug_id, target_symbol, target_weight, is_direct, action_type."""
    dt = drug_targets.merge(universe[["drug_id"]], on="drug_id", how="inner")
    dt = dt.drop_duplicates(["drug_id", "target_symbol"])
    dt["target_weight"] = 1.0
    dt["is_direct"] = True
    cols = ["drug_id", "target_symbol", "target_weight", "is_direct", "action_type"]

    if not cfg.get("network", "enabled", default=False):
        return dt[cols]

    client = StringClient(
        source=cfg.get("network", "source", default="file"),
        edges_path=cfg.path("string_edges") if cfg.get("network", "source") == "file" else None,
        min_confidence=cfg.get("network", "min_confidence", default=0.7),
        max_partners=cfg.get("network", "max_partners", default=10),
        species=cfg.get("network", "string_species", default=9606),
        cache_dir=cfg.root / cfg.get("network", "cache_dir", default=".cache/string"),
    )
    nbr_w = float(cfg.get("network", "neighbor_weight", default=0.5))
    partners = client.neighbors(sorted(dt["target_symbol"].unique().tolist()))
    if partners.empty:
        return dt[cols]

    exp = dt[["drug_id", "target_symbol", "action_type"]].merge(partners, on="target_symbol", how="inner")
    exp["target_symbol"] = exp["partner_symbol"]
    exp["target_weight"] = nbr_w * exp["score"]
    exp["is_direct"] = False
    # action_type is only meaningful for *direct* targets, so neighbour edges carry no
    # reliable pharmacological direction.
    exp["action_type"] = ""
    exp = exp[cols]

    both = pd.concat([dt[cols], exp])
    both = (both.sort_values("target_weight", ascending=False)
                .drop_duplicates(["drug_id", "target_symbol"]))  # keep the stronger (direct) edge
    return both.reset_index(drop=True)


# ----------------------------------------------------------------------
# Step 3a - drug x target x disease edges (shared substrate)
# ----------------------------------------------------------------------
def disease_edges(drug_targets: pd.DataFrame, target_disease: pd.DataFrame,
                  cfg: Config) -> pd.DataFrame:
    min_assoc = float(cfg.get("propagation", "min_assoc", default=0.1))
    td = target_disease[target_disease["assoc_score"] >= min_assoc]
    e = drug_targets.merge(td, on="target_symbol", how="inner")
    if e.empty:
        return e.assign(contrib=[])
    e["contrib"] = e["target_weight"].clip(0, 1) * e["assoc_score"].clip(0, 1)
    return e


# ----------------------------------------------------------------------
# Step 3b - aggregate edges -> mechanistic support per (drug, disease)
# ----------------------------------------------------------------------
def propagate_disease(edges: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    agg = cfg.get("propagation", "aggregate", default="noisy_or")
    keep = int(cfg.get("propagation", "top_targets_kept", default=5))
    cols = ["drug_id", "efo_id", "disease_name", "mechanistic_support", "n_targets",
            "lead_target", "evidence_targets"]
    if edges.empty:
        return pd.DataFrame(columns=cols)

    def _aggregate(g: pd.DataFrame) -> pd.Series:
        c = g["contrib"].values
        if agg == "max":
            m = float(c.max())
        elif agg == "sum":
            m = float(min(c.sum(), 1.0))
        else:
            m = float(1.0 - np.prod(1.0 - np.clip(c, 0, 0.999)))
        g_sorted = g.sort_values("contrib", ascending=False)
        top = g_sorted.head(keep).apply(lambda r: f"{r['target_symbol']}({r['contrib']:.2f})", axis=1)
        return pd.Series({
            "disease_name": g["disease_name"].iloc[0],
            "mechanistic_support": m,
            "n_targets": int(g["target_symbol"].nunique()),
            "lead_target": g_sorted.iloc[0]["target_symbol"],
            "evidence_targets": ", ".join(top),
        })

    return (edges.groupby(["drug_id", "efo_id"], group_keys=True)
                 .apply(_aggregate, include_groups=False).reset_index())[cols]


# ----------------------------------------------------------------------
# Step 3c - directionality model
# ----------------------------------------------------------------------
_ACTION_DOWN = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR",
                "NEGATIVE MODULATOR", "INVERSE AGONIST", "DISRUPTING AGENT", "DEGRADER",
                "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "SUPPRESSOR"}
_ACTION_UP = {"AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR",
              "PARTIAL AGONIST", "OPENER", "STABILISER", "STABILIZER"}


def action_to_dir(action_type: str) -> int:
    """+1 if the drug raises target activity, -1 if it lowers it, 0 if unknown."""
    a = str(action_type).upper().strip()
    if a in _ACTION_DOWN:
        return -1
    if a in _ACTION_UP:
        return 1
    return 0


def add_direction(edges: pd.DataFrame, target_direction: pd.DataFrame,
                  cfg: Config) -> pd.DataFrame:
    """Per (drug, disease): direction_factor, direction_status, direction_alignment.

    alignment(edge) = drug_dir(action_type) * therapeutic_direction
        +1 = drug pushes the target the therapeutically useful way (aligned)
        -1 = drug pushes it the wrong way (opposed)
    Only DIRECT targets with a known action and known disease direction contribute.
    The per-disease alignment is the contribution-weighted mean of edge alignments.
    """
    aligned = float(cfg.get("direction", "aligned_factor", default=1.0))
    opposed = float(cfg.get("direction", "opposed_factor", default=0.15))
    unknown = float(cfg.get("direction", "unknown_factor", default=0.6))
    thr = float(cfg.get("direction", "align_threshold", default=0.34))

    out_cols = ["drug_id", "efo_id", "direction_factor", "direction_status", "direction_alignment"]
    if edges.empty:
        return pd.DataFrame(columns=out_cols)

    e = edges[edges["is_direct"]].copy()
    e["drug_dir"] = e["action_type"].map(action_to_dir)
    e = e.merge(target_direction[["target_symbol", "efo_id", "therapeutic_direction"]],
                on=["target_symbol", "efo_id"], how="left")
    e["therapeutic_direction"] = e["therapeutic_direction"].fillna(0).astype(int)
    e["informative"] = (e["drug_dir"] != 0) & (e["therapeutic_direction"] != 0)
    e["alignment"] = e["drug_dir"] * e["therapeutic_direction"]

    def _factor(g: pd.DataFrame) -> pd.Series:
        gi = g[g["informative"]]
        if gi.empty or gi["contrib"].sum() == 0:
            return pd.Series({"direction_factor": unknown, "direction_status": "unknown",
                              "direction_alignment": np.nan})
        al = float((gi["contrib"] * gi["alignment"]).sum() / gi["contrib"].sum())
        factor = opposed + (aligned - opposed) * (al + 1.0) / 2.0
        status = "aligned" if al >= thr else ("opposed" if al <= -thr else "mixed")
        return pd.Series({"direction_factor": round(factor, 4),
                          "direction_status": status, "direction_alignment": round(al, 3)})

    res = (e.groupby(["drug_id", "efo_id"], group_keys=True)
            .apply(_factor, include_groups=False).reset_index())
    return res[out_cols]


# ----------------------------------------------------------------------
# Step 3d - tissue-expression filter
# ----------------------------------------------------------------------
def add_tissue(edges: pd.DataFrame, target_expression: pd.DataFrame,
               disease_tissue: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Per (drug, disease): tissue_factor, tissue_status, tissue_evidence.

    Asks 'is a driving target actually expressed where the disease manifests?'.
    For each supporting target we take its max expression across the disease's
    relevant tissues; the hypothesis takes the best (max) over its targets.
        expressed (>= min_expression) -> tissue_factor = 1.0
        low (0 < score < min)         -> scaled down toward absent_factor
        absent (target measured, 0)   -> absent_factor
        unknown (no expression data)  -> unknown_factor (mild discount, not a veto)
    """
    min_expr = float(cfg.get("tissue", "min_expression", default=0.25))
    absent = float(cfg.get("tissue", "absent_factor", default=0.3))
    unknown = float(cfg.get("tissue", "unknown_factor", default=0.7))
    out_cols = ["drug_id", "efo_id", "tissue_factor", "tissue_status", "tissue_evidence"]
    if edges.empty:
        return pd.DataFrame(columns=out_cols)

    measured = set(target_expression["target_symbol"])
    mapped = set(disease_tissue["efo_id"])

    # tissue_expr(target, disease) = max over relevant tissues of relevance * expression
    te = disease_tissue.merge(target_expression, on="tissue", how="inner")
    te["te"] = te["relevance"].clip(0, 1) * te["expression"].clip(0, 1)
    te_max = (te.sort_values("te", ascending=False)
                .drop_duplicates(["target_symbol", "efo_id"])
                [["target_symbol", "efo_id", "te", "tissue"]])

    e = edges[["drug_id", "efo_id", "target_symbol"]].merge(
        te_max, on=["target_symbol", "efo_id"], how="left")
    e["informative"] = e["efo_id"].isin(mapped) & e["target_symbol"].isin(measured)
    e["te_edge"] = np.where(e["informative"], e["te"].fillna(0.0), np.nan)

    def _factor(g: pd.DataFrame) -> pd.Series:
        gi = g.dropna(subset=["te_edge"])
        if gi.empty:
            return pd.Series({"tissue_factor": unknown, "tissue_status": "unknown",
                              "tissue_evidence": ""})
        best = gi.loc[gi["te_edge"].idxmax()]
        score_ = float(best["te_edge"])
        tissue = best["tissue"] if isinstance(best["tissue"], str) else ""
        if score_ >= min_expr:
            return pd.Series({"tissue_factor": 1.0, "tissue_status": "expressed", "tissue_evidence": tissue})
        if score_ > 0:
            f = round(absent + (1.0 - absent) * score_ / min_expr, 4)
            return pd.Series({"tissue_factor": f, "tissue_status": "low", "tissue_evidence": tissue})
        return pd.Series({"tissue_factor": absent, "tissue_status": "absent", "tissue_evidence": ""})

    res = (e.groupby(["drug_id", "efo_id"], group_keys=True)
            .apply(_factor, include_groups=False).reset_index())
    return res[out_cols]


# ----------------------------------------------------------------------
# Step 3e - phylogenetics (orthologous-gene model-organism evidence)
# ----------------------------------------------------------------------
def add_phylo_evidence(edges: pd.DataFrame, phylo: pd.DataFrame,
                       cfg: Config) -> pd.DataFrame:
    """Per (drug, disease): phylo_factor, phylo_score, phylo_n_models, phylo_sources.

    Asymmetric -- present evidence BOOSTS the opportunity score, absence
    does NOT penalise (factor 1.0 in both 'no measurement' and 'measured
    no signal'). Rationale: many target-disease relationships in the
    literature ARE evidenced through orthologous knockouts in mouse/fly/
    fish, but the absence of such a study is not evidence against the
    target -- the model organism may not have been studied yet.

    Aggregation across the drug's targets is MAX -- a single target with
    a strong cross-species phenotype is sufficient to credit the drug.

        factor = 1 + boost_factor * phylo_score   (clip to [1, 1 + boost_factor])
        boost_factor = 0   -> no effect (inspection-only via the column)
        boost_factor = 0.5 -> up to 1.5x on opportunity (default)
    """
    boost = float(cfg.get("phylogenetics", "boost_factor", default=0.5))
    out_cols = ["drug_id", "efo_id", "phylo_factor", "phylo_score",
                "phylo_n_models", "phylo_sources"]
    if edges.empty or phylo.empty:
        return pd.DataFrame(columns=out_cols)

    e = edges[["drug_id", "efo_id", "target_symbol"]].merge(
        phylo[["target_symbol", "efo_id", "phylo_score", "n_models", "sources"]],
        on=["target_symbol", "efo_id"], how="inner")
    if e.empty:
        return pd.DataFrame(columns=out_cols)

    # Per (drug, disease): take the target with the strongest phylo_score
    # (one target with strong model-organism evidence is enough to boost).
    idx = e.groupby(["drug_id", "efo_id"])["phylo_score"].idxmax()
    best = e.loc[idx, ["drug_id", "efo_id", "phylo_score", "n_models", "sources"]].copy()
    best = best.rename(columns={"n_models": "phylo_n_models",
                                "sources": "phylo_sources"})
    best["phylo_factor"] = (1.0 + boost * best["phylo_score"].clip(0, 1)).round(4)
    return best[out_cols].reset_index(drop=True)


# ----------------------------------------------------------------------
# Step 4a - close ChEMBL indication-coverage gaps via target-class roll-up
# ----------------------------------------------------------------------
def expand_indications_by_target_class(
    indications: pd.DataFrame, drug_targets: pd.DataFrame, cfg: Config,
) -> pd.DataFrame:
    """Augment indications by pooling them across drugs in the same target class.

    ChEMBL annotates indications on individual drug entries, but many drugs
    that ARE in the same pharmacological class get patchy coverage:
      * formulation variants (INSULIN ASPART vs INSULIN ASPART PROTAMINE
        RECOMBINANT) -- same exact targets, different ChEMBL rows
      * selective vs broad-spectrum members of a class (silodosin is
        ADRA1A-selective; doxazosin is pan-ADRA1) -- silodosin's class is
        a SUBSET of doxazosin's, so silodosin should inherit doxazosin's
        hypertension indication (we want to suppress 'use alpha-1 blocker
        for HTN' as a 'novel' hypothesis -- it's textbook).

    Roll-up rule: drug A inherits drug B's indications iff A's frozenset
    of (direct target_symbol, action_type) pairs is a SUBSET of B's.
    Equality is the trivial case; strict-subset is the selective-inherits-
    from-broad case. STRING neighbours are excluded so they don't inflate
    sharing. A kinase inhibitor with a unique target set still doesn't
    inherit from a different kinase inhibitor (neither is a subset of
    the other).
    """
    if not bool(cfg.get("novelty", "expand_by_target_class", default=True)):
        return indications
    if indications.empty or drug_targets.empty:
        return indications
    dt = drug_targets.copy()
    if "is_direct" in dt.columns:
        dt = dt[dt["is_direct"]]
    if dt.empty:
        return indications
    dt["pair"] = list(zip(dt["target_symbol"].fillna(""),
                          dt["action_type"].fillna("")))
    drug_class = (dt.groupby("drug_id")["pair"]
                    .apply(lambda s: frozenset(s)).rename("class_set").reset_index())

    # Indications grouped by the contributor drug's class set (set of efo_ids).
    ind_with_class = indications.merge(drug_class, on="drug_id", how="inner")
    class_to_efos: dict[frozenset, set[str]] = {}
    for cls, efo in zip(ind_with_class["class_set"], ind_with_class["efo_id"]):
        class_to_efos.setdefault(cls, set()).add(efo)

    # For each unique drug class C, pool indications from every class C' with C ⊆ C'.
    contributor_classes = list(class_to_efos.keys())
    pooled_by_class: dict[frozenset, set[str]] = {}
    for cls in set(drug_class["class_set"]):
        pooled: set[str] = set()
        for c2 in contributor_classes:
            if cls <= c2:                      # equality handled here too
                pooled.update(class_to_efos[c2])
        if pooled:
            pooled_by_class[cls] = pooled

    # Distribute back to drugs.
    rolled_rows = [
        {"drug_id": row.drug_id, "efo_id": efo, "indication_name": ""}
        for row in drug_class.itertuples(index=False)
        for efo in pooled_by_class.get(row.class_set, ())
    ]
    rolled = (pd.DataFrame(rolled_rows) if rolled_rows
              else pd.DataFrame(columns=["drug_id", "efo_id", "indication_name"]))
    combined = pd.concat([indications, rolled], ignore_index=True).drop_duplicates(
        subset=["drug_id", "efo_id"], keep="first").reset_index(drop=True)
    return combined


# ----------------------------------------------------------------------
# Step 4 - novelty
# ----------------------------------------------------------------------
def known_expanded(indications: pd.DataFrame, ontology: pd.DataFrame,
                   cfg: Config) -> pd.DataFrame:
    """Drug-disease pairs that are KNOWN (or ontology-near a known indication).

    Returns drug_id, efo_id, novelty, novelty_status for those pairs only; any
    pair absent from this table is treated as fully novel (novelty 1.0).
    """
    radius = int(cfg.get("novelty", "ontology_radius", default=1))
    soft = bool(cfg.get("novelty", "soft", default=True))
    onto = DiseaseOntology(ontology)

    best: dict[tuple[str, str], tuple[float, str]] = {}
    for drug_id, grp in indications.groupby("drug_id"):
        known = set(grp["efo_id"])
        nbhd = onto.neighborhood(known, radius) if radius > 0 else {e: 0 for e in known}
        for efo, hop in nbhd.items():
            if hop == 0:
                nov, st = 0.0, "known_exact"
            else:
                nov = round(hop / (radius + 1), 3) if soft else 0.0
                st = f"known_related_h{hop}"
            key = (drug_id, efo)
            if key not in best or nov < best[key][0]:
                best[key] = (nov, st)
    if not best:
        return pd.DataFrame(columns=["drug_id", "efo_id", "novelty", "novelty_status"])
    rows = [(d, e, nov, st) for (d, e), (nov, st) in best.items()]
    return pd.DataFrame(rows, columns=["drug_id", "efo_id", "novelty", "novelty_status"])


def add_novelty(hypotheses: pd.DataFrame, known_exp: pd.DataFrame) -> pd.DataFrame:
    out = hypotheses.merge(known_exp, on=["drug_id", "efo_id"], how="left")
    out["novelty"] = out["novelty"].fillna(1.0)
    out["novelty_status"] = out["novelty_status"].fillna("novel")
    return out


# ----------------------------------------------------------------------
# Step 5 - opportunity score + ranking
# ----------------------------------------------------------------------
def score(hypotheses: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    w_m = float(cfg.get("scoring", "w_mech", default=1.0))
    w_n = float(cfg.get("scoring", "w_novelty", default=1.0))
    penalty = bool(cfg.get("scoring", "promiscuity_penalty", default=True))
    min_opp = float(cfg.get("scoring", "min_opportunity", default=0.0))
    top_n = int(cfg.get("scoring", "top_n", default=100000))

    h = hypotheses.copy()
    if "direction_factor" not in h.columns:
        h["direction_factor"] = float(cfg.get("direction", "default_factor", default=1.0))
    h["direction_factor"] = h["direction_factor"].fillna(
        float(cfg.get("direction", "unknown_factor", default=0.6)))
    if "tissue_factor" not in h.columns:
        h["tissue_factor"] = 1.0
    h["tissue_factor"] = h["tissue_factor"].fillna(
        float(cfg.get("tissue", "unknown_factor", default=0.7)))
    # phylo_factor is BOOST-ONLY: missing or absent -> 1.0 (no penalty).
    if "phylo_factor" not in h.columns:
        h["phylo_factor"] = 1.0
    h["phylo_factor"] = h["phylo_factor"].fillna(1.0)

    opp = ((h["mechanistic_support"].clip(0, 1) ** w_m)
           * (h["novelty"].clip(0, 1) ** w_n)
           * h["direction_factor"]
           * h["tissue_factor"]
           * h["phylo_factor"])
    if penalty:
        basis = h["n_drug_targets"] if "n_drug_targets" in h.columns else h["n_targets"]
        opp = opp / np.sqrt(basis.clip(lower=1))
    h["opportunity"] = opp.round(5)
    h = h[h["opportunity"] >= min_opp]
    h = h.sort_values("opportunity", ascending=False).head(top_n).reset_index(drop=True)
    h.insert(0, "rank", h.index + 1)
    return h


# ----------------------------------------------------------------------
# Step 5b - two-pass literature / patent / trial enrichment
# ----------------------------------------------------------------------
def add_literature_pass(ranked: pd.DataFrame, lit_client, cfg: Config) -> pd.DataFrame:
    """Query the literature client for the top-K hypotheses, damp by prior, re-rank.

    Two-pass design avoids hitting external APIs for every (drug, disease) pair:
      1. The earlier scoring pass produced `ranked` sorted by opportunity.
      2. We take `literature.top_n` rows, dedup on (lead_target, disease_name),
         query the client (cache covers most calls on re-runs), merge back.
      3. Damp opportunity by  (1 - w_investigation * investigation_prior)  and
         re-sort + re-rank.

    Rows outside the top-K keep their original opportunity / rank and get NaN
    literature columns. Rows inside the top-K but with a missing
    (lead_target, disease_name) (rare; nothing was looked up) also keep their
    original opportunity but get NaN literature columns.
    """
    top_n = int(cfg.get("literature", "top_n", default=5000))
    w_inv = float(cfg.get("scoring", "w_investigation", default=0.5))
    count_cols = ["pubmed_count", "europepmc_count", "trial_count", "patent_count"]
    out = ranked.copy()
    if out.empty or top_n <= 0:
        for c in count_cols:
            out[c] = pd.array([pd.NA] * len(out), dtype="Int64")
        out["investigation_prior"] = float("nan")
        return out

    head = out.head(top_n)
    uniq_cols = ["lead_target", "efo_id", "disease_name"]
    for c in ("target_synonyms", "disease_synonyms"):
        if c in head.columns:
            uniq_cols.append(c)
    uniq = (head[uniq_cols]
            .rename(columns={"lead_target": "target_symbol"})
            .dropna(subset=["target_symbol", "disease_name"])
            .drop_duplicates(subset=["target_symbol", "efo_id", "disease_name"]))
    prior = lit_client.score_pairs(uniq).rename(columns={"target_symbol": "lead_target"})
    # Drop the synonym helper columns from the merge result -- they were query
    # aids; the only thing we keep is the count columns + investigation_prior.
    drop_helpers = [c for c in ("target_synonyms", "disease_synonyms") if c in prior.columns]
    if drop_helpers:
        prior = prior.drop(columns=drop_helpers)
    out = out.merge(prior, on=["lead_target", "efo_id", "disease_name"], how="left")
    if w_inv > 0:
        damp = (1.0 - w_inv * out["investigation_prior"].fillna(0.0).clip(0, 1)).clip(0, 1)
        # Only damp rows that actually received a literature lookup -- a missing
        # prior means "we didn't query it", not "we queried and found nothing".
        looked_up = out["investigation_prior"].notna()
        out.loc[looked_up, "opportunity"] = (
            out.loc[looked_up, "opportunity"] * damp.loc[looked_up]
        ).round(5)
        out = out.sort_values("opportunity", ascending=False).reset_index(drop=True)
        if "rank" in out.columns:
            out["rank"] = out.index + 1
    return out


# ----------------------------------------------------------------------
# Step 5f - pathway-level mechanistic evidence (Reactome via OT)
# ----------------------------------------------------------------------
def add_pathway_evidence(ranked: pd.DataFrame, drug_targets: pd.DataFrame,
                          target_pathways: pd.DataFrame,
                          target_disease: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Boost opportunity when drug's targets share pathways with disease's strong targets.

    Captures the indirect mechanism case: drug hits target A; target A is in
    the same pathway as disease-driving target B that the drug doesn't hit
    directly. Reactome pathway co-membership is treated as weaker evidence
    than direct target overlap (which mechanistic_support already credits),
    so the default boost is smaller.

    Restricted to SPECIFIC pathways (<= max_pathway_size genes) to avoid
    generic categories matching everything. Restricted to STRONGLY-
    associated disease targets (assoc_score >= min_assoc) to avoid
    pulling in disease associations driven by noise.

    Asymmetric boost like phylogenetics: present overlap multiplies
    opportunity by (1 + boost * pathway_score); absence applies no penalty
    (factor 1.0).

        pathway_score = log1p(n_overlap) / log1p(saturation)
        pathway_factor = 1 + boost * pathway_score   (default boost = 0.3)
    """
    if not bool(cfg.get("pathway", "enabled", default=False)):
        ranked = ranked.copy()
        ranked["pathway_factor"] = 1.0
        ranked["pathway_score"] = 0.0
        ranked["n_pathway_overlap"] = 0
        return ranked
    boost = float(cfg.get("pathway", "boost_factor", default=0.3))
    max_size = int(cfg.get("pathway", "max_pathway_size", default=80))
    min_assoc = float(cfg.get("pathway", "min_assoc", default=0.3))
    saturation = int(cfg.get("pathway", "saturation", default=5))
    indirect_only = bool(cfg.get("pathway", "indirect_only", default=True))

    out = ranked.copy()
    if target_pathways.empty or out.empty:
        out["pathway_factor"] = 1.0
        out["pathway_score"] = 0.0
        out["n_pathway_overlap"] = 0
        return out

    # Filter pathways to specific ones (small enough to be meaningful)
    sizes = target_pathways.groupby("pathway_id")["target_symbol"].nunique()
    keep = set(sizes[sizes <= max_size].index)
    tp = target_pathways[target_pathways["pathway_id"].isin(keep)]

    # disease_pathways: pathway memberships of the disease's STRONG targets,
    # keeping target identity so we can filter direct hits later.
    td_strong = target_disease[target_disease["assoc_score"] >= min_assoc]
    disease_paths = (td_strong[["efo_id", "target_symbol"]]
                       .merge(tp[["target_symbol", "pathway_id"]], on="target_symbol")
                       .drop_duplicates(["efo_id", "target_symbol", "pathway_id"])
                       .rename(columns={"target_symbol": "disease_target"}))

    # drug_pathways: pathway memberships of the drug's DIRECT targets.
    dt = drug_targets
    if "is_direct" in dt.columns:
        dt = dt[dt["is_direct"]]
    drug_paths = (dt[["drug_id", "target_symbol"]]
                    .merge(tp[["target_symbol", "pathway_id"]], on="target_symbol")
                    .drop_duplicates(["drug_id", "target_symbol", "pathway_id"])
                    .rename(columns={"target_symbol": "drug_target"}))

    # Bridge: for each (drug, disease, pathway), find drug_target and
    # disease_target sharing that pathway.
    keys = out[["drug_id", "efo_id"]].drop_duplicates()
    bridges = (keys.merge(drug_paths, on="drug_id")
                   .merge(disease_paths, on=["efo_id", "pathway_id"]))
    if indirect_only:
        # The pathway boost should reward cases the direct-target signal MISSES.
        # Drop bridges where the drug hits exactly the disease's strong target
        # in this pathway -- that contribution is already in mechanistic_support.
        bridges = bridges[bridges["drug_target"] != bridges["disease_target"]]
    overlap = (bridges.drop_duplicates(["drug_id", "efo_id", "pathway_id"])
                       .groupby(["drug_id", "efo_id"]).size()
                       .rename("n_pathway_overlap").reset_index())
    overlap["pathway_score"] = overlap["n_pathway_overlap"].apply(
        lambda c: round(min(np.log1p(c) / np.log1p(saturation), 1.0), 4))
    overlap["pathway_factor"] = (1.0 + boost * overlap["pathway_score"]).round(4)

    out = out.merge(overlap, on=["drug_id", "efo_id"], how="left")
    out["pathway_factor"] = out["pathway_factor"].fillna(1.0)
    out["pathway_score"] = out["pathway_score"].fillna(0.0)
    out["n_pathway_overlap"] = out["n_pathway_overlap"].fillna(0).astype(int)
    # Apply the multiplicative boost to opportunity and re-rank
    out["opportunity"] = (out["opportunity"] * out["pathway_factor"]).round(5)
    out = out.sort_values("opportunity", ascending=False).reset_index(drop=True)
    if "rank" in out.columns:
        out["rank"] = out.index + 1
    return out


# ----------------------------------------------------------------------
# Step 5h - severity flag: receptor-LoF + agonist drug = futile
# ----------------------------------------------------------------------
_SEVERITY_KEYWORDS = (
    "deficiency", "complete absence", "complete loss", "null mutation",
    "null variant", "donohue", "rabson-mendenhall",
)


def flag_severity_concern(ranked: pd.DataFrame, drug_targets: pd.DataFrame,
                          cfg: Config) -> pd.DataFrame:
    """Flag (drug, disease) pairs where a `disease_gene_match` row also has
    severity language ('deficiency', etc.) AND the drug is an AGONIST on the
    named gene -- the classic 'give the agonist to the broken receptor' trap.

    INSULIN (INSR AGONIST) -> 'hyperinsulinism due to INSR deficiency' (futile)
    PEGFILGRASTIM (CSF3R AGONIST) -> '...due to CSF3R deficiency' (futile)

    But: TZDs -> 'PPARG-related familial partial lipodystrophy' (no
    'deficiency' keyword -> not flagged -> kept as a legit hypothesis).
    Calcimimetics -> 'familial hypocalciuric hypercalcemia 1' (same).

    Optionally damps opportunity by severity.damp_factor (default 0.7 = a
    71% knockdown). Inspection-only when damp_factor == 0.
    """
    if not bool(cfg.get("severity", "enabled", default=True)):
        out = ranked.copy()
        out["severity_concern"] = ""
        return out

    damp = float(cfg.get("severity", "damp_factor", default=0.7))
    out = ranked.copy()
    out["severity_concern"] = ""
    if "disease_gene_match" not in out.columns or "disease_name" not in out.columns:
        return out

    # Per (drug, target) action_type lookup over DIRECT targets only.
    dt = drug_targets
    if "is_direct" in dt.columns:
        dt = dt[dt["is_direct"]]
    actions: dict[tuple[str, str], str] = {
        (row.drug_id, row.target_symbol): (row.action_type or "").upper().strip()
        for row in dt[["drug_id", "target_symbol", "action_type"]].itertuples(index=False)
    }
    agonist_actions = _ACTION_UP | {"AGONIST", "ACTIVATOR"}

    flags = []
    for _, row in out.iterrows():
        match = row.get("disease_gene_match") or ""
        if not match:
            flags.append("")
            continue
        disease_name = (row.get("disease_name") or "").lower()
        if not any(kw in disease_name for kw in _SEVERITY_KEYWORDS):
            flags.append("")
            continue
        triggered = False
        for gene in match.split(","):
            action = actions.get((row.drug_id, gene), "")
            if action in agonist_actions:
                triggered = True
                break
        flags.append("severe_loF_agonist" if triggered else "")
    out["severity_concern"] = flags

    if damp > 0 and (out["severity_concern"] != "").any():
        mask = out["severity_concern"] == "severe_loF_agonist"
        out.loc[mask, "opportunity"] = (
            out.loc[mask, "opportunity"] * (1.0 - damp)
        ).round(5)
        out = out.sort_values("opportunity", ascending=False).reset_index(drop=True)
        if "rank" in out.columns:
            out["rank"] = out.index + 1
    return out


# ----------------------------------------------------------------------
# Step 5e - active-substance grouping (collapse formulation variants)
# ----------------------------------------------------------------------
def collapse_to_substances(ranked: pd.DataFrame, substance_map: pd.DataFrame,
                            cfg: Config) -> pd.DataFrame:
    """Collapse (drug_id, efo_id) rows into (substance_chembl_id, efo_id).

    Without this, the output is dominated by formulation variants of the
    same active ingredient -- 11 INSULIN rows for hyperinsulinism, 4 TZD
    rows for FPLD3, etc. -- all making the same hypothesis with identical
    scores. We collapse by ChEMBL's active-ingredient parent and keep the
    highest-opportunity row per substance-disease, while preserving the
    variant names as a comma-separated `variant_names` column.

    Mode flags via config.substance_grouping:
      enabled: false -> pass through unchanged (no new columns)
      enabled, collapse=false -> add substance columns but no row collapse
      enabled, collapse=true  -> collapse + re-rank (default)
    """
    if not bool(cfg.get("substance_grouping", "enabled", default=True)):
        return ranked
    out = ranked.copy()
    if substance_map is None or substance_map.empty:
        out["substance_chembl_id"] = out["drug_id"]
        out["substance_name"] = out.get("drug_name", out["drug_id"])
        out["n_variants"] = 1
        out["variant_names"] = out.get("drug_name", out["drug_id"])
        return out
    out = out.merge(substance_map, on="drug_id", how="left")
    out["substance_chembl_id"] = out["substance_chembl_id"].fillna(out["drug_id"])
    if "drug_name" in out.columns:
        out["substance_name"] = out["substance_name"].fillna(out["drug_name"])
    else:
        out["substance_name"] = out["substance_name"].fillna(out["substance_chembl_id"])

    if not bool(cfg.get("substance_grouping", "collapse", default=True)):
        out["n_variants"] = 1
        out["variant_names"] = out.get("drug_name", out["substance_name"])
        return out

    # Collapse: highest-opportunity row per (substance, disease) is the
    # representative; variant_names lists every formulation that mapped
    # to this substance for this disease.
    out = out.sort_values("opportunity", ascending=False, kind="stable")
    variant_col = "drug_name" if "drug_name" in out.columns else "drug_id"
    variants = (out.groupby(["substance_chembl_id", "efo_id"])[variant_col]
                  .apply(lambda s: list(dict.fromkeys(s)))
                  .rename("variants_list").reset_index())
    variants["n_variants"] = variants["variants_list"].apply(len)
    variants["variant_names"] = variants["variants_list"].apply(lambda lst: ", ".join(lst))
    best = out.drop_duplicates(subset=["substance_chembl_id", "efo_id"], keep="first")
    best = best.merge(variants[["substance_chembl_id", "efo_id",
                                "n_variants", "variant_names"]],
                      on=["substance_chembl_id", "efo_id"], how="left")
    # Replace the displayed drug_name with the canonical substance_name
    if "drug_name" in best.columns:
        best["drug_name"] = best["substance_name"]
    best = best.sort_values("opportunity", ascending=False).reset_index(drop=True)
    if "rank" in best.columns:
        best["rank"] = best.index + 1
    return best


# ----------------------------------------------------------------------
# Step 5g - KOL finder (one US + one EU key opinion leader per hypothesis)
# ----------------------------------------------------------------------
def add_kol_pass(ranked: pd.DataFrame, kol_client, cfg: Config) -> pd.DataFrame:
    """Look up a US + EU KOL for the top-N (lead_target, disease) pairs.

    Inspection-only: never changes opportunity / rank. Uses the same
    pattern as the literature pass -- dedup on (lead_target, efo_id,
    disease_name) so substance variants share lookups.
    """
    top_n = int(cfg.get("kol", "top_n", default=100))
    out_cols = [f"{r}_kol_{c}" for r in ("us", "eu")
                for c in ("name", "institution", "email", "h_index", "n_pubs")]
    out = ranked.copy()
    if out.empty or top_n <= 0:
        for c in out_cols:
            out[c] = "" if not c.endswith(("h_index", "n_pubs")) else pd.NA
        return out

    head = out.head(top_n)
    uniq_cols = ["lead_target", "efo_id", "disease_name"]
    if "disease_synonyms" in head.columns:
        uniq_cols.append("disease_synonyms")
    uniq = (head[uniq_cols]
            .rename(columns={"lead_target": "target_symbol"})
            .dropna(subset=["target_symbol", "disease_name"])
            .drop_duplicates(subset=["target_symbol", "efo_id", "disease_name"]))
    kols = kol_client.find_kols(uniq).rename(columns={"target_symbol": "lead_target"})
    # Strip synonyms column from the merge result so we don't duplicate it on ranked
    if "disease_synonyms" in kols.columns and "disease_synonyms" in out.columns:
        kols = kols.drop(columns=["disease_synonyms"])
    out = out.merge(kols, on=["lead_target", "efo_id", "disease_name"], how="left")
    for c in out_cols:
        if c.endswith(("h_index", "n_pubs")):
            out[c] = out[c].astype("Int64") if c in out.columns else pd.NA
        else:
            out[c] = out[c].fillna("") if c in out.columns else ""
    return out


# ----------------------------------------------------------------------
# Step 5d - flag (drug, disease) pairs where the drug targets a gene named in the disease
# ----------------------------------------------------------------------
_GENE_TOKEN_RE = re.compile(r"\b([A-Z][A-Z0-9]{1,6})\b")
_RECEPTOR_SUFFIX_RE = re.compile(r"\b([A-Z][A-Z0-9]{1,6})\s+receptor\b")


def _disease_gene_tokens(name: str, gene_universe: set[str]) -> set[str]:
    """Extract gene-like uppercase tokens from a disease name, restricted to
    symbols that actually exist in the gene universe. Also derives 'XR' from
    'X receptor' (TSH receptor -> TSHR, EPO receptor -> EPOR) to catch the
    common naming convention where the ligand is named but the receptor is
    the actual target.
    """
    if not isinstance(name, str) or not name:
        return set()
    tokens = set(_GENE_TOKEN_RE.findall(name))
    for m in _RECEPTOR_SUFFIX_RE.finditer(name):
        tokens.add(m.group(1) + "R")
    return tokens & gene_universe


def flag_disease_gene_concordance(hypotheses: pd.DataFrame,
                                  drug_targets: pd.DataFrame,
                                  gene_info: pd.DataFrame) -> pd.DataFrame:
    """Per (drug, disease): which of the drug's direct targets are named in the disease.

    Surfaces the whole class of 'drug targets the disease gene directly' hits:
      * INSULIN -> 'hyperinsulinism due to INSR deficiency'         (severe LoF -- agonist futile)
      * THYROTROPIN -> 'hypothyroidism due to TSH receptor mutations' (severe LoF -- agonist futile)
      * PIOGLITAZONE -> 'PPARG-related familial partial lipodystrophy' (partial LoF -- legit)
      * CINACALCET -> 'familial hypocalciuric hypercalcemia 1'        (CASR; need full name to match)

    Inspection-only: emits a `disease_gene_match` column (comma-separated
    symbols, empty when none). Whether the match means the drug is wrong
    (severe receptor LoF -> agonists can't help) or right (partial LoF
    rescuable by agonists) depends on disease-specific severity -- defer
    to human judgment rather than auto-damp.
    """
    out_cols = ["drug_id", "efo_id", "disease_gene_match"]
    if hypotheses.empty:
        return pd.DataFrame(columns=out_cols)

    gene_universe: set[str] = set()
    if gene_info is not None and not gene_info.empty and "symbol" in gene_info.columns:
        gene_universe = set(gene_info["symbol"].dropna().astype(str))
    if drug_targets is not None and not drug_targets.empty:
        gene_universe |= set(drug_targets["target_symbol"].dropna().astype(str))
    if not gene_universe:
        return hypotheses[["drug_id", "efo_id"]].assign(disease_gene_match="")

    # Per drug -> set of direct target symbols.
    dt = drug_targets.copy()
    if "is_direct" in dt.columns:
        dt = dt[dt["is_direct"].fillna(True)]
    drug_to_targets = (dt.dropna(subset=["target_symbol"])
                         .groupby("drug_id")["target_symbol"]
                         .apply(lambda s: set(s.astype(str))).to_dict())

    # Per disease name -> set of gene-symbol tokens (cache by name; many drugs share).
    unique_names = hypotheses[["efo_id", "disease_name"]].drop_duplicates()
    disease_tokens = {
        efo_id: _disease_gene_tokens(name, gene_universe)
        for efo_id, name in zip(unique_names["efo_id"], unique_names["disease_name"])
    }

    matches = []
    for drug_id, efo_id in zip(hypotheses["drug_id"], hypotheses["efo_id"]):
        toks = disease_tokens.get(efo_id, set())
        if not toks:
            matches.append("")
            continue
        targets = drug_to_targets.get(drug_id, set())
        hit = sorted(toks & targets)
        matches.append(",".join(hit))
    return hypotheses[["drug_id", "efo_id"]].assign(disease_gene_match=matches)


# ----------------------------------------------------------------------
# Step 5c - US market-size enrichment (inspection-only)
# ----------------------------------------------------------------------
def add_market_size(ranked: pd.DataFrame, market_client, cfg: Config) -> pd.DataFrame:
    """Attach US market data to the top-N hypotheses without changing ranking.

    Single-pass over unique EFO IDs in the top-N (small fan-out: ~thousands
    of distinct diseases even at full scale). The orphan-drug flag is
    derived from `us_patients` against the US Orphan Drug Act threshold
    (default 200,000); rows whose patient count is unknown get NA, NOT
    False -- we don't want to silently say "not orphan" when we don't know.
    """
    top_n = int(cfg.get("market", "top_n", default=10_000))
    orphan_threshold = int(cfg.get("market", "rare_disease_us_threshold", default=200_000))
    out = ranked.copy()
    if out.empty or top_n <= 0:
        out["us_patients"] = pd.array([pd.NA] * len(out), dtype="Int64")
        out["us_prevalence_per_100k"] = float("nan")
        out["market_source"] = ""
        out["as_of"] = ""
        out["is_orphan"] = pd.array([pd.NA] * len(out), dtype="boolean")
        return out

    head = out.head(top_n)
    uniq_ids = head["efo_id"].dropna().drop_duplicates().tolist()
    market = market_client.lookup(uniq_ids)
    out = out.merge(market, on="efo_id", how="left")
    out["us_patients"] = out["us_patients"].astype("Int64")
    out["market_source"] = out["market_source"].fillna("")
    out["as_of"] = out["as_of"].fillna("")
    out["is_orphan"] = pd.array(
        [(x < orphan_threshold) if pd.notna(x) else pd.NA for x in out["us_patients"]],
        dtype="boolean")
    return out
