"""DuckDB backend for the full-scale join.

The big relation is drug x target x disease. Instead of loading the entire
Open Targets association table into pandas, DuckDB scans it from disk (CSV or
Parquet) and performs the join, the noisy-OR aggregation, the directionality
aggregation, the tissue-expression gate, the novelty join and the final scoring
in a single SQL pass.

All small, pre-computed inputs (drug-target edges, directionality evidence,
expression / tissue maps, known-indication closure, drug breadth) are registered
as in-memory relations. The output schema matches the pandas engine.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

from ..config import Config
from ..steps import _ACTION_DOWN, _ACTION_UP


def _sql_in(values) -> str:
    return ", ".join("'" + v.replace("'", "''") + "'" for v in sorted(values))


def _td_scan(path: Path) -> str:
    p = str(path).replace("'", "''")
    if Path(path).suffix.lower() in {".parquet", ".pq"}:
        return f"parquet_scan('{p}')"
    return f"read_csv_auto('{p}', header=true, quote='\"', escape='\"', sample_size=-1)"


def run_duckdb(cfg: Config, prep: dict, log=lambda m: None) -> pd.DataFrame:
    import duckdb

    min_assoc = float(cfg.get("propagation", "min_assoc", default=0.1))
    agg = cfg.get("propagation", "aggregate", default="noisy_or")
    w_m = float(cfg.get("scoring", "w_mech", default=1.0))
    w_n = float(cfg.get("scoring", "w_novelty", default=1.0))
    penalty = bool(cfg.get("scoring", "promiscuity_penalty", default=True))
    min_opp = float(cfg.get("scoring", "min_opportunity", default=0.0))
    top_n = int(cfg.get("scoring", "top_n", default=100000))

    direction_on = bool(prep["direction_on"])
    aligned = float(cfg.get("direction", "aligned_factor", default=1.0))
    opposed = float(cfg.get("direction", "opposed_factor", default=0.15))
    d_unknown = float(cfg.get("direction", "unknown_factor", default=0.6))
    thr = float(cfg.get("direction", "align_threshold", default=0.34))
    default_factor = float(cfg.get("direction", "default_factor", default=1.0))

    tissue_on = bool(prep["tissue_on"])
    min_expr = float(cfg.get("tissue", "min_expression", default=0.25))
    absent = float(cfg.get("tissue", "absent_factor", default=0.3))
    t_unknown = float(cfg.get("tissue", "unknown_factor", default=0.7))

    phylo_on = bool(prep.get("phylo_on", False))
    phylo_boost = float(cfg.get("phylogenetics", "boost_factor", default=0.5))

    con = duckdb.connect()
    # Hint DuckDB to be modest with memory/threads -- the disease_tissue x
    # target_expression cross-product is hot enough to thrash without bounds.
    # Fewer threads = smaller multiplicative temp footprint per stage.
    con.execute("PRAGMA memory_limit='8GB'")
    con.execute("PRAGMA threads=2")
    con.execute("PRAGMA preserve_insertion_order=false")
    con.register("dt", prep["drug_targets"])
    con.register("tdir", prep["target_direction"])
    con.register("texpr", prep["target_expression"])
    con.register("dtissue", prep["disease_tissue"])
    con.register("phylo_ev", prep.get("phylo_evidence", pd.DataFrame(
        columns=["target_symbol", "efo_id", "phylo_score", "n_models", "sources"])))
    con.register("known_exp", prep["known_exp"])
    con.register("breadth", prep["breadth"])

    support_expr = {
        "max": "MAX(contrib)",
        "sum": "LEAST(SUM(contrib), 1.0)",
    }.get(agg, "1 - EXP(SUM(LN(1 - LEAST(contrib, 0.999))))")
    down, up = _sql_in(_ACTION_DOWN), _sql_in(_ACTION_UP)
    pen_expr = "SQRT(GREATEST(COALESCE(b.n_drug_targets, 1), 1))" if penalty else "1.0"

    # ---- direction ----
    if direction_on:
        dir_cte = f"""
        , dir AS (
          SELECT drug_id, efo_id,
                 SUM(contrib*drug_dir*ther_dir) FILTER (WHERE drug_dir<>0 AND ther_dir<>0) AS num,
                 SUM(contrib)                    FILTER (WHERE drug_dir<>0 AND ther_dir<>0) AS den
          FROM (
            SELECT e.drug_id, e.efo_id, e.contrib,
                   CASE WHEN UPPER(TRIM(e.action_type)) IN ({down}) THEN -1
                        WHEN UPPER(TRIM(e.action_type)) IN ({up})   THEN  1 ELSE 0 END AS drug_dir,
                   COALESCE(td.therapeutic_direction, 0) AS ther_dir
            FROM edges e LEFT JOIN tdir td
              ON e.target_symbol = td.target_symbol AND e.efo_id = td.efo_id
            WHERE e.is_direct
          ) GROUP BY 1, 2
        )"""
        dir_factor = (f"ROUND(CASE WHEN d.den IS NULL OR d.den = 0 THEN {d_unknown} "
                      f"ELSE {opposed} + ({aligned}-{opposed})*((d.num/d.den)+1)/2 END, 4)")
        dir_status = (f"CASE WHEN d.den IS NULL OR d.den = 0 THEN 'unknown' "
                      f"WHEN d.num/d.den >= {thr} THEN 'aligned' "
                      f"WHEN d.num/d.den <= -{thr} THEN 'opposed' ELSE 'mixed' END")
        dir_join = "LEFT JOIN dir d ON s.drug_id = d.drug_id AND s.efo_id = d.efo_id"
    else:
        dir_cte, dir_factor, dir_status, dir_join = "", f"{default_factor}", "NULL", ""

    # ---- tissue ----
    if tissue_on:
        tissue_cte = f"""
        , te_map AS (
          SELECT efo_id, target_symbol, MAX(rel_expr) AS te, arg_max(tissue, rel_expr) AS best_tissue
          FROM (
            SELECT dd.efo_id, tt.target_symbol, dd.tissue,
                   LEAST(dd.relevance,1.0)*LEAST(tt.expression,1.0) AS rel_expr
            FROM dtissue dd JOIN texpr tt USING (tissue)
          ) GROUP BY 1, 2
        ),
        tissue AS (
          SELECT drug_id, efo_id, MAX(te_edge) AS tscore, COUNT(te_edge) AS n_inf,
                 arg_max(best_tissue, te_edge) AS best_tissue
          FROM (
            SELECT e.drug_id, e.efo_id,
                   CASE WHEN e.efo_id IN (SELECT efo_id FROM dtissue)
                             AND e.target_symbol IN (SELECT target_symbol FROM texpr)
                        THEN COALESCE(tm.te, 0.0) ELSE NULL END AS te_edge,
                   tm.best_tissue
            FROM edges e LEFT JOIN te_map tm
              ON e.target_symbol = tm.target_symbol AND e.efo_id = tm.efo_id
          ) GROUP BY 1, 2
        )"""
        tis_factor = (f"ROUND(CASE WHEN t.n_inf IS NULL OR t.n_inf = 0 THEN {t_unknown} "
                      f"WHEN t.tscore >= {min_expr} THEN 1.0 "
                      f"WHEN t.tscore > 0 THEN {absent} + (1-{absent})*t.tscore/{min_expr} "
                      f"ELSE {absent} END, 4)")
        tis_status = (f"CASE WHEN t.n_inf IS NULL OR t.n_inf = 0 THEN 'unknown' "
                      f"WHEN t.tscore >= {min_expr} THEN 'expressed' "
                      f"WHEN t.tscore > 0 THEN 'low' ELSE 'absent' END")
        tis_evidence = "CASE WHEN t.tscore > 0 THEN COALESCE(t.best_tissue, '') ELSE '' END"
        tis_join = "LEFT JOIN tissue t ON s.drug_id = t.drug_id AND s.efo_id = t.efo_id"
    else:
        tissue_cte, tis_factor, tis_status, tis_evidence, tis_join = "", "1.0", "NULL", "''", ""

    # ---- phylogenetics (orthologous-gene model-organism evidence) ----
    if phylo_on:
        phylo_cte = """
        , phylo AS (
          SELECT e.drug_id, e.efo_id,
                 MAX(p.phylo_score) AS phylo_score,
                 arg_max(p.n_models, p.phylo_score) AS phylo_n_models,
                 arg_max(p.sources, p.phylo_score) AS phylo_sources
          FROM edges e JOIN phylo_ev p
            ON e.target_symbol = p.target_symbol AND e.efo_id = p.efo_id
          GROUP BY 1, 2
        )"""
        # Boost-only: COALESCE to 1.0 when there's no phylo row (no penalty).
        phy_factor = (f"ROUND(COALESCE(1.0 + {phylo_boost} * "
                      f"LEAST(GREATEST(ph.phylo_score, 0), 1), 1.0), 4)")
        phy_score = "COALESCE(ph.phylo_score, 0.0)"
        phy_models = "COALESCE(ph.phylo_n_models, 0)"
        phy_sources = "COALESCE(ph.phylo_sources, '')"
        phy_join = "LEFT JOIN phylo ph ON s.drug_id = ph.drug_id AND s.efo_id = ph.efo_id"
    else:
        phylo_cte, phy_factor = "", "1.0"
        phy_score, phy_models, phy_sources, phy_join = "NULL", "NULL", "NULL", ""

    query = f"""
    WITH td AS (
        SELECT CAST(target_symbol AS VARCHAR) target_symbol,
               CAST(efo_id AS VARCHAR) efo_id,
               CAST(disease_name AS VARCHAR) disease_name,
               CAST(assoc_score AS DOUBLE) assoc_score
        FROM {_td_scan(cfg.path('target_disease'))}
        WHERE CAST(assoc_score AS DOUBLE) >= {min_assoc}
    ),
    edges AS (
        SELECT dt.drug_id, dt.target_symbol, dt.action_type, dt.is_direct,
               td.efo_id, td.disease_name,
               LEAST(dt.target_weight, 1.0) * LEAST(td.assoc_score, 1.0) AS contrib
        FROM dt JOIN td USING (target_symbol)
    ),
    support AS (
        SELECT drug_id, efo_id, ANY_VALUE(disease_name) AS disease_name,
               {support_expr} AS mechanistic_support,
               COUNT(DISTINCT target_symbol) AS n_targets,
               arg_max(target_symbol, contrib) AS lead_target,
               string_agg(target_symbol || '(' || ROUND(contrib, 2) || ')', ', '
                          ORDER BY contrib DESC) AS evidence_targets
        FROM edges GROUP BY drug_id, efo_id
    )
    {dir_cte}
    {tissue_cte}
    {phylo_cte}
    , scored AS (
        SELECT s.drug_id, s.efo_id, s.disease_name, s.lead_target,
               s.mechanistic_support,
               COALESCE(k.novelty, 1.0) AS novelty,
               COALESCE(k.novelty_status, 'novel') AS novelty_status,
               {dir_factor} AS direction_factor,
               {dir_status} AS direction_status,
               {tis_factor} AS tissue_factor,
               {tis_status} AS tissue_status,
               {tis_evidence} AS tissue_evidence,
               {phy_factor} AS phylo_factor,
               {phy_score} AS phylo_score,
               {phy_models} AS phylo_n_models,
               {phy_sources} AS phylo_sources,
               s.n_targets,
               COALESCE(b.n_drug_targets, 1) AS n_drug_targets,
               ROUND( POWER(LEAST(GREATEST(s.mechanistic_support,0),1), {w_m})
                    * POWER(LEAST(GREATEST(COALESCE(k.novelty,1.0),0),1), {w_n})
                    * ({dir_factor}) * ({tis_factor}) * ({phy_factor}) / {pen_expr}, 5) AS opportunity,
               s.evidence_targets
        FROM support s
        LEFT JOIN known_exp k ON s.drug_id = k.drug_id AND s.efo_id = k.efo_id
        LEFT JOIN breadth   b ON s.drug_id = b.drug_id
        {dir_join}
        {tis_join}
        {phy_join}
    )
    SELECT * FROM scored
    WHERE opportunity >= {min_opp}
    ORDER BY opportunity DESC
    LIMIT {top_n}
    """
    df = con.execute(query).fetch_df()
    con.close()
    df.insert(0, "rank", range(1, len(df) + 1))
    log(f"duckdb produced {len(df)} ranked hypotheses")
    return df
