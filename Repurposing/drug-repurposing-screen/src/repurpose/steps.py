"""Pipeline steps. Each is a pure function over pandas DataFrames.

The expensive object is `disease_edges` (drug x target x disease). Both the
mechanistic-support aggregation and the directionality model are computed from
it, so the two engines (pandas / duckdb) stay consistent.
"""
from __future__ import annotations
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

    opp = ((h["mechanistic_support"].clip(0, 1) ** w_m)
           * (h["novelty"].clip(0, 1) ** w_n)
           * h["direction_factor"]
           * h["tissue_factor"])
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
    lit_cols = ["pubmed_count", "europepmc_count", "trial_count", "investigation_prior"]
    out = ranked.copy()
    if out.empty or top_n <= 0:
        for c in lit_cols[:3]:
            out[c] = pd.array([pd.NA] * len(out), dtype="Int64")
        out["investigation_prior"] = float("nan")
        return out

    head = out.head(top_n)
    uniq = (head[["lead_target", "efo_id", "disease_name"]]
            .rename(columns={"lead_target": "target_symbol"})
            .dropna(subset=["target_symbol", "disease_name"])
            .drop_duplicates())
    prior = lit_client.score_pairs(uniq).rename(columns={"target_symbol": "lead_target"})
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
