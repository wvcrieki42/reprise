"""Adapters that turn raw public bulk dumps into the canonical CSV/Parquet tables.

These are the concrete extractors for `mode: full`. Run them once to materialise the
five canonical tables, then point config.yaml at the outputs.

Data sources (all free / open):
  * ChEMBL    : approved drugs, mechanism-of-action targets, drug indications
                https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/  (chembl_*_sqlite.tar.gz)
  * Open Targets associations : target -> disease scores (parquet)
                https://platform.opentargets.org/downloads  (associationByOverallDirect)
  * EFO       : disease ontology parents (efo.obo / owl) for novelty radius
  * DrugBank  : optional cross-check of approved status / targets (XML, licensed)

Each function returns a pandas DataFrame in canonical schema. They are intentionally
self-contained and dependency-light (sqlite3, pandas) so they run on a laptop.
"""
from __future__ import annotations
import sqlite3
import pandas as pd


# ----------------------------------------------------------------------
# ChEMBL (SQLite)  ->  drugs, drug_targets, drug_indications
# ----------------------------------------------------------------------
def chembl_drugs(sqlite_path: str) -> pd.DataFrame:
    """Approved/launched drugs with region flags inferred from regulatory fields."""
    sql = """
    SELECT md.chembl_id              AS drug_id,
           md.pref_name              AS drug_name,
           md.max_phase              AS max_phase,
           CASE WHEN md.usan_stem IS NOT NULL OR md.first_approval IS NOT NULL
                THEN 1 ELSE 0 END    AS approved_us,
           0                         AS approved_eu,
           CASE WHEN md.molecule_type = 'Small molecule' THEN 'small_molecule'
                WHEN md.molecule_type = 'Antibody'       THEN 'antibody'
                WHEN md.molecule_type = 'Protein'        THEN 'protein'
                WHEN md.molecule_type = 'Oligonucleotide' THEN 'oligonucleotide'
                ELSE 'other' END     AS modality
    FROM molecule_dictionary md
    WHERE md.max_phase >= 4;
    """
    con = sqlite3.connect(sqlite_path)
    try:
        df = pd.read_sql(sql, con)
    finally:
        con.close()
    # NOTE: ChEMBL alone is US-centric for approval. Merge EMA/Orange Book lists
    # (see scripts/download_data.sh) to set approved_eu accurately.
    return df


def chembl_drug_targets(sqlite_path: str) -> pd.DataFrame:
    """Curated mechanism-of-action targets (the drug's *known* MoA)."""
    sql = """
    SELECT md.chembl_id      AS drug_id,
           cs.component_synonym AS target_symbol,
           dm.action_type    AS action_type,
           dm.mechanism_of_action AS mechanism_of_action
    FROM drug_mechanism dm
    JOIN molecule_dictionary md ON md.molregno = dm.molregno
    JOIN target_components tc   ON tc.tid = dm.tid
    JOIN component_synonyms cs  ON cs.component_id = tc.component_id
                               AND cs.syn_type = 'GENE_SYMBOL'
    WHERE md.max_phase >= 4;
    """
    con = sqlite3.connect(sqlite_path)
    try:
        df = pd.read_sql(sql, con).dropna(subset=["target_symbol"]).drop_duplicates()
    finally:
        con.close()
    return df


def chembl_substance_map(sqlite_path: str) -> pd.DataFrame:
    """drug_id -> substance_chembl_id, substance_name.

    Rolls each drug up to its ChEMBL active-ingredient parent via
    molecule_hierarchy.parent_molregno. Catches the small-molecule salt
    case cleanly (PIOGLITAZONE HYDROCHLORIDE -> PIOGLITAZONE, KETAMINE
    HYDROCHLORIDE -> KETAMINE). For biologics where ChEMBL stores each
    formulation as its own root entry (INSULIN ASPART vs INSULIN ASPART
    PROTAMINE RECOMBINANT both self-parented), the substance_chembl_id
    equals the drug_id -- a name-based fallback would be needed to
    collapse those further; not yet wired.
    """
    sql = """
    SELECT md.chembl_id              AS drug_id,
           COALESCE(pmd.chembl_id, md.chembl_id) AS substance_chembl_id,
           COALESCE(pmd.pref_name,  md.pref_name) AS substance_name
    FROM molecule_dictionary md
    LEFT JOIN molecule_hierarchy mh ON mh.molregno = md.molregno
    LEFT JOIN molecule_dictionary pmd ON pmd.molregno = mh.parent_molregno
    WHERE md.chembl_id IS NOT NULL
    """
    con = sqlite3.connect(sqlite_path)
    try:
        df = pd.read_sql(sql, con).drop_duplicates(subset=["drug_id"])
    finally:
        con.close()
    return df


def chembl_drug_indications(sqlite_path: str) -> pd.DataFrame:
    """Known indications (subtracted from the screen to leave novel hypotheses)."""
    sql = """
    SELECT md.chembl_id   AS drug_id,
           di.efo_id      AS efo_id,
           di.efo_term    AS indication_name
    FROM drug_indication di
    JOIN molecule_dictionary md ON md.molregno = di.molregno
    WHERE di.efo_id IS NOT NULL;
    """
    con = sqlite3.connect(sqlite_path)
    try:
        df = pd.read_sql(sql, con).drop_duplicates()
    finally:
        con.close()
    # ChEMBL exports ontology IDs in colon form for non-EFO sources
    # (MONDO:0005148, HP:0000863, etc.). The rest of the pipeline uses
    # the underscore form (MONDO_0005148, HP_0000863, EFO_0001360). Without
    # this normalisation, MONDO-coded indications NEVER match OT's
    # target_disease EFO IDs and approved drugs (SGLT2 inhibitors for T2D,
    # CFTR modulators for cystic fibrosis, etc.) keep getting flagged
    # "novel" because the novelty subtraction silently fails.
    df["efo_id"] = df["efo_id"].str.replace(":", "_", regex=False)
    return df


# ----------------------------------------------------------------------
# Open Targets associations (parquet)  ->  target_disease
# ----------------------------------------------------------------------
def opentargets_target_disease(parquet_dir: str, gene_map_csv: str,
                               disease_map_csv: str | None = None,
                               min_score: float = 0.1) -> pd.DataFrame:
    """target -> disease association scores.

    parquet_dir : Open Targets `associationByOverallDirect` parquet directory.
    gene_map_csv: two-column map ensembl_id,target_symbol (from OT `targets` dataset).
    """
    assoc = pd.read_parquet(parquet_dir)
    assoc = assoc.rename(columns={
        "targetId": "ensembl_id", "diseaseId": "efo_id", "score": "assoc_score"})
    assoc = assoc[assoc["assoc_score"] >= min_score]
    genes = pd.read_csv(gene_map_csv)  # ensembl_id, target_symbol
    out = assoc.merge(genes, on="ensembl_id", how="inner")
    if disease_map_csv:
        dm = pd.read_csv(disease_map_csv)
        if {"efo_id", "disease_name"}.issubset(dm.columns):
            dm = dm[["efo_id", "disease_name"]].dropna(subset=["efo_id"]).drop_duplicates(subset=["efo_id"])
            out = out.merge(dm, on="efo_id", how="left")
        else:
            out["disease_name"] = ""
    else:
        out["disease_name"] = ""
    out["disease_name"] = out["disease_name"].fillna("")
    return out[["target_symbol", "efo_id", "disease_name", "assoc_score"]]


# ----------------------------------------------------------------------
# EFO ontology  ->  disease_ontology (efo_id, name, parent_efo_id)
# ----------------------------------------------------------------------
# OT datasources whose variantEffect / directionOnTrait fields actually encode
# THERAPEUTIC direction in humans. Mouse-KO phenotypes (IMPC) and somatic-cancer
# sources are intentionally EXCLUDED here -- they are great for `add_phylo_evidence`
# but not for inferring agonist-vs-antagonist therapeutic intent. Without this
# filter, IMPC alone (~1.1M rows, the largest source by far) drowns out the
# curated human-genetics signals and flips therapeutic direction for monogenic
# disorders where the human variant mechanism (e.g. truncating GoF EPOR) differs
# from what a complete mouse knockout looks like.
DEFAULT_DIRECTION_SOURCES: tuple[str, ...] = (
    "ot_genetics_portal",   # GWAS effect-direction
    "gene_burden",          # rare-variant burden tests
    "eva",                  # ClinVar germline interpretations
    "gene2phenotype",       # curated gene-disease
    "clingen",              # ClinGen gene-disease validity
    "genomics_england",     # PanelApp
    "orphanet",             # rare-disease curation
)


def opentargets_target_direction(evidence_parquet_dir: str, gene_map_csv: str,
                                 sources: tuple[str, ...] = DEFAULT_DIRECTION_SOURCES) -> pd.DataFrame:
    """Derive therapeutic_direction in {-1,+1} from Open Targets genetic evidence.

    Uses the variant effect (loss/gain of function) and its direction on the trait
    (risk/protective) to infer which way a drug should move the target:

        increasing target activity is RISK       -> inhibit  (therapeutic_direction = -1)
        increasing target activity is PROTECTIVE -> activate (therapeutic_direction = +1)

    where 'increasing activity' = GoF-risk / LoF-protective etc. Conflicting
    genes are resolved by majority vote across the kept sources. Only the OT
    sources listed in `sources` (default DEFAULT_DIRECTION_SOURCES) are read --
    we explicitly exclude IMPC and other non-therapeutic sources here.
    """
    from pathlib import Path as _Path
    files: list[str] = []
    base = _Path(evidence_parquet_dir)
    for src in sources:
        files.extend(str(p) for p in sorted(base.glob(f"sourceId={src}/*.parquet")))
    if not files:
        return pd.DataFrame(columns=["target_symbol", "efo_id",
                                     "therapeutic_direction", "evidence"])
    ev = pd.read_parquet(files, columns=[
        "targetId", "diseaseId", "variantEffect", "directionOnTrait"])
    ev = ev.rename(columns={"targetId": "ensembl_id", "diseaseId": "efo_id"})
    ev = ev.dropna(subset=["variantEffect", "directionOnTrait"])

    def _therapeutic_dir(row) -> int:
        eff = str(row["variantEffect"]).upper()        # LOF / GOF
        trait = str(row["directionOnTrait"]).upper()    # RISK / PROTECTIVE
        # sign of how *increasing activity* affects disease:
        #   GoF-risk  or LoF-protective -> increasing activity is harmful  -> inhibit (-1)
        #   GoF-protective or LoF-risk  -> increasing activity is helpful  -> activate (+1)
        inc_harmful = (("GOF" in eff and "RISK" in trait) or ("LOF" in eff and "PROTECT" in trait))
        inc_helpful = (("GOF" in eff and "PROTECT" in trait) or ("LOF" in eff and "RISK" in trait))
        return -1 if inc_harmful else (1 if inc_helpful else 0)

    ev["d"] = ev.apply(_therapeutic_dir, axis=1)
    ev = ev[ev["d"] != 0]
    vote = (ev.groupby(["ensembl_id", "efo_id"])["d"].sum().reset_index())
    vote["therapeutic_direction"] = vote["d"].apply(lambda x: 1 if x > 0 else (-1 if x < 0 else 0))
    vote = vote[vote["therapeutic_direction"] != 0]
    genes = pd.read_csv(gene_map_csv)  # ensembl_id, target_symbol
    out = vote.merge(genes, on="ensembl_id", how="inner")
    out["evidence"] = "Open Targets genetic direction-of-effect"
    return out[["target_symbol", "efo_id", "therapeutic_direction", "evidence"]]


def opentargets_phylo_evidence(evidence_parquet_dir: str, gene_map_csv: str,
                               sources: tuple[str, ...] = ("impc",),
                               min_score: float = 0.0) -> pd.DataFrame:
    """Per (target_symbol, efo_id): model-organism (orthologous gene) evidence.

    Reads OT evidence rows where `datatypeId == "animal_model"` from the
    requested sources (default: IMPC mouse knockouts; PhenoDigm provides
    the cross-species phenotype matching that produces the score). Each
    row is a piece of evidence -- a specific mouse model exhibiting a
    phenotype that maps to the human disease.

    Aggregation:
      * per (target_symbol, efo_id): take MAX score across rows
        (strongest evidence wins; many weak hits != one strong cross-
        species match)
      * also report n_models (number of evidence rows kept) and the set
        of contributing sources, for downstream provenance.

    Returns columns: target_symbol, efo_id, phylo_score, n_models, sources.
    """
    ev = pd.read_parquet(evidence_parquet_dir, columns=[
        "datasourceId", "datatypeId", "targetId", "diseaseId", "score"])
    ev = ev[(ev["datatypeId"] == "animal_model") & (ev["datasourceId"].isin(sources))]
    ev = ev[ev["score"] >= min_score]
    if ev.empty:
        return pd.DataFrame(columns=["target_symbol", "efo_id", "phylo_score",
                                     "n_models", "sources"])
    ev = ev.rename(columns={"targetId": "ensembl_id", "diseaseId": "efo_id"})
    agg = (ev.groupby(["ensembl_id", "efo_id"])
             .agg(phylo_score=("score", "max"),
                  n_models=("score", "size"),
                  sources=("datasourceId", lambda s: ",".join(sorted(set(s)))))
             .reset_index())
    genes = pd.read_csv(gene_map_csv)  # ensembl_id, target_symbol
    out = agg.merge(genes, on="ensembl_id", how="inner")
    return out[["target_symbol", "efo_id", "phylo_score", "n_models", "sources"]]


def opentargets_gene_info(targets_parquet_dir: str) -> pd.DataFrame:
    """symbol -> full target name, from the Open Targets `targets` dataset."""
    t = pd.read_parquet(targets_parquet_dir)
    t = t.rename(columns={"approvedSymbol": "symbol", "approvedName": "gene_name"})
    return t[["symbol", "gene_name"]].dropna().drop_duplicates()


def opentargets_target_expression(baseline_parquet_dir: str, gene_map_csv: str,
                                  rna_field: str = "rna") -> pd.DataFrame:
    """target x tissue baseline expression, normalised to [0,1].

    Parses the Open Targets `baselineExpression` dataset (per-target tissue records)
    and rescales RNA level to 0..1 by the per-target maximum. `gene_map_csv` maps
    ensembl_id -> target_symbol. Adapt the column names to the OT release in use.
    """
    be = pd.read_parquet(baseline_parquet_dir)  # columns: id (ensembl), tissues (list of structs)
    rows = []
    for ensembl, tissues in zip(be["id"], be["tissues"]):
        if tissues is None or (isinstance(tissues, float) and pd.isna(tissues)):
            tissue_entries = []
        elif isinstance(tissues, (list, tuple)):
            tissue_entries = tissues
        else:
            try:
                tissue_entries = list(tissues)
            except TypeError:
                tissue_entries = []
        for ts in tissue_entries:
            lvl = ts.get(rna_field) if isinstance(ts, dict) else None
            if isinstance(lvl, dict):
                lvl = lvl.get("value")
            if lvl is not None:
                rows.append((ensembl, ts.get("label", ts.get("efo_code", "")), float(lvl)))
    df = pd.DataFrame(rows, columns=["ensembl_id", "tissue", "level"])
    df["expression"] = df.groupby("ensembl_id")["level"].transform(
        lambda s: s / s.max() if s.max() else 0.0)
    genes = pd.read_csv(gene_map_csv)  # ensembl_id, target_symbol
    out = df.merge(genes, on="ensembl_id", how="inner")
    return out[["target_symbol", "tissue", "expression"]]


def efo_ontology(obo_path: str) -> pd.DataFrame:
    """Parse an EFO .obo into (efo_id, disease_name, parent_efo_id) edges."""
    rows, cur, name = [], None, None
    with open(obo_path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                cur, name = None, None
            elif line.startswith("id: "):
                cur = line[4:].replace(":", "_")
            elif line.startswith("name: "):
                name = line[6:]
            elif line.startswith("is_a: ") and cur:
                parent = line[6:].split("!")[0].strip().replace(":", "_")
                rows.append((cur, name or "", parent))
    df = pd.DataFrame(rows, columns=["efo_id", "disease_name", "parent_efo_id"])
    # ensure every node appears at least once even with no parent
    return df.drop_duplicates()
