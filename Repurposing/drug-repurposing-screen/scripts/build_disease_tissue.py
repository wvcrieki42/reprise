"""Bootstrap data/full/disease_tissue.csv from OT diseases parquet.

For each disease, emit one row per (disease, tissue, relevance) based on
its OT therapeutic areas mapped to anatomical tissues via the curated
TA_TO_TISSUE dict. On top of that, a smaller hand-curated
DISEASE_NAME_OVERRIDES dict adds disease-specific tissues (cancers are
the main case -- a cancer's primary site can't be inferred from the
"cancer or benign tumor" therapeutic area alone).

Tissue names use the OT baseline-expression vocabulary EXACTLY (e.g.
"heart left ventricle", "skeletal muscle tissue", "prostate gland") so
the disease_tissue.csv joins cleanly to target_expression.csv on the
`tissue` column.

Existing rows in disease_tissue.csv are PRESERVED on (efo_id, tissue)
collision -- manual edits win on re-run.

Behaviour around overly-broad therapeutic areas:
  * Cancer (MONDO_0045024) intentionally has NO general tissue. We
    don't want every cancer entry tagged with a default site. The
    DISEASE_NAME_OVERRIDES catches the major cancers.
  * Genetic/familial/congenital (OTAR_0000018), phenotype, measurement,
    biological_process, injury/poisoning -- skipped, no anatomy signal.

Run after scripts/download_data.sh and scripts/build_full_tables.py.
"""
from __future__ import annotations
import glob
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DISEASES_DIR = ROOT / "data" / "full" / "diseases"
OUT_CSV = ROOT / "data" / "full" / "disease_tissue.csv"


# Therapeutic area EFO/MONDO id  ->  list of (tissue, relevance)
# Tissue names MUST match the OT baseline-expression vocabulary in
# data/full/target_expression.csv exactly.
TA_TO_TISSUE: dict[str, list[tuple[str, float]]] = {
    # cardiovascular disease
    "EFO_0000319": [("heart", 1.0), ("heart left ventricle", 0.9),
                    ("heart muscle", 0.9), ("aorta", 0.7), ("coronary artery", 0.7)],
    # hematologic disease
    "EFO_0005803": [("blood", 1.0), ("bone marrow", 0.9)],
    # immune system disease
    "EFO_0000540": [("blood", 0.8), ("spleen", 0.8), ("thymus", 0.7),
                    ("bone marrow", 0.6)],
    # nervous system disease
    "EFO_0000618": [("brain", 1.0), ("cerebral cortex", 0.8), ("cerebellum", 0.6)],
    # psychiatric disorder
    "MONDO_0002025": [("brain", 1.0), ("cerebral cortex", 0.9)],
    # respiratory or thoracic disease
    "OTAR_0000010": [("lung", 1.0), ("bronchus", 0.8)],
    # endocrine system disease
    "EFO_0010282": [("thyroid gland", 0.7), ("adrenal gland", 0.7),
                    ("pancreas", 0.5)],
    "EFO_0001379": [("thyroid gland", 0.7), ("adrenal gland", 0.7),
                    ("pancreas", 0.5)],
    # gastrointestinal disease
    "EFO_0010285": [("colon", 0.9), ("small intestine", 0.9), ("stomach", 0.8),
                    ("esophagus", 0.6), ("duodenum", 0.7)],
    # musculoskeletal / connective tissue disease
    "OTAR_0000006": [("skeletal muscle tissue", 0.9)],
    # integumentary (skin)
    "MONDO_0024458": [("eye", 1.0)],  # MONDO disorder of visual system
    # urinary system disease
    "EFO_0009690": [("kidney", 1.0), ("cortex of kidney", 0.9)],
    # reproductive system or breast disease
    "OTAR_0000017": [("breast", 0.7), ("ovary", 0.7), ("uterus", 0.7),
                     ("prostate gland", 0.7)],
    # pancreas disease
    "EFO_0009605": [("pancreas", 1.0)],
    # disorder of ear
    "MONDO_0021205": [("inner ear", 1.0)] if False else [],   # 'ear' is not in OT expression vocab; skip
    # nutritional / metabolic
    "OTAR_0000020": [("liver", 0.6), ("adipose tissue", 0.6)],
}

# Skin doesn't have its own TA; integumentary-system-disease falls under nervous TA
# in OT's classification. We add it via name overrides instead.

# Disease-name-specific overrides (lowercase exact match against OT name).
# Cancer primary-site tagging + a few common-condition refinements.
DISEASE_NAME_OVERRIDES: dict[str, list[tuple[str, float]]] = {
    # cancer primary sites
    "breast cancer": [("breast", 1.0)],
    "prostate cancer": [("prostate gland", 1.0)],
    "colorectal cancer": [("colon", 1.0)],
    "lung cancer": [("lung", 1.0)],
    "ovarian cancer": [("ovary", 1.0)],
    "pancreatic carcinoma": [("pancreas", 1.0)],
    "pancreatic adenocarcinoma": [("pancreas", 1.0)],
    "thyroid cancer": [("thyroid gland", 1.0)],
    "uterine cancer": [("uterus", 1.0)],
    "kidney cancer": [("kidney", 1.0)],
    "renal cell carcinoma": [("kidney", 1.0), ("cortex of kidney", 0.9)],
    "urinary bladder cancer": [("urinary bladder", 1.0)] if False else [],   # not in OT vocab
    "gastric cancer": [("stomach", 1.0)],
    "esophageal cancer": [("esophagus", 1.0)],
    "hepatocellular carcinoma": [("liver", 1.0)],
    "cervical cancer": [("uterus", 0.7)],
    "head and neck squamous cell carcinoma": [("tongue", 0.8), ("esophagus", 0.5)],
    "gastrointestinal stromal tumor": [("stomach", 0.8), ("small intestine", 0.8)],
    "melanoma": [("skin of body", 1.0)],
    # leukemias / lymphomas already covered by hematologic TA, but reinforce
    "chronic lymphocytic leukemia": [("blood", 1.0), ("bone marrow", 0.9)],
    "chronic myelogenous leukemia": [("blood", 1.0), ("bone marrow", 0.9)],
    "acute myeloid leukemia": [("blood", 1.0), ("bone marrow", 0.9)],
    "acute lymphoblastic leukemia": [("blood", 1.0), ("bone marrow", 0.9)],
    "multiple myeloma": [("bone marrow", 1.0), ("blood", 0.9)],
    "glioblastoma multiforme": [("brain", 1.0), ("cerebral cortex", 0.8)],
    "glioma": [("brain", 1.0), ("cerebral cortex", 0.8)],
    # common conditions where the TA mapping needs reinforcement
    "idiopathic pulmonary fibrosis": [("lung", 1.0)],
    "asthma": [("lung", 1.0), ("bronchus", 0.9)],
    "chronic obstructive pulmonary disease": [("lung", 1.0), ("bronchus", 0.9)],
    "atopic eczema": [("skin of body", 1.0)],
    "psoriasis": [("skin of body", 1.0)],
    "rosacea": [("skin of body", 1.0)],
    "vitiligo": [("skin of body", 1.0)],
    "alopecia areata": [("skin of body", 0.9)],
    "hidradenitis suppurativa": [("skin of body", 1.0)],
    "ulcerative colitis": [("colon", 1.0)],
    "Crohn's disease": [("small intestine", 1.0), ("colon", 0.8)],
    "inflammatory bowel disease": [("colon", 1.0), ("small intestine", 1.0)],
    "type 2 diabetes mellitus": [("pancreas", 1.0), ("liver", 0.7),
                                 ("adipose tissue", 0.6)],
    "type 1 diabetes mellitus": [("pancreas", 1.0)],
    "alzheimer disease": [("brain", 1.0), ("cerebral cortex", 1.0)],
    "parkinson disease": [("brain", 1.0)],
    "amyotrophic lateral sclerosis": [("spinal cord", 1.0), ("brain", 0.7)] if False else [("brain", 0.9)],  # spinal cord not in vocab as plain term
    "huntington disease": [("brain", 1.0)],
    "multiple sclerosis": [("brain", 1.0), ("cerebral cortex", 0.7)],
    "lewy body dementia": [("brain", 1.0)],
    "epilepsy": [("brain", 1.0)],
    "migraine disorder": [("brain", 1.0)],
    "non-alcoholic fatty liver disease": [("liver", 1.0)],
    "non-alcoholic steatohepatitis": [("liver", 1.0)],
    "chronic kidney disease": [("kidney", 1.0)],
    "iga glomerulonephritis": [("kidney", 1.0), ("cortex of kidney", 1.0)],
    "autosomal dominant polycystic kidney disease": [("kidney", 1.0)],
    "rheumatoid arthritis": [("skeletal muscle tissue", 0.5)],   # synovium not in vocab; placeholder
    "ankylosing spondylitis": [("skeletal muscle tissue", 0.5)],
    "psoriatic arthritis": [("skeletal muscle tissue", 0.5), ("skin of body", 0.9)],
    "age-related macular degeneration": [("retina", 1.0)] if False else [("eye", 1.0)],
    "glaucoma": [("eye", 1.0)],
    "diabetic retinopathy": [("eye", 1.0)],
    "retinitis pigmentosa": [("eye", 1.0)],
    "atrial fibrillation": [("heart", 1.0), ("atrium auricular region", 1.0)],
    "heart failure": [("heart", 1.0), ("heart left ventricle", 1.0)],
    "coronary artery disease": [("heart", 1.0), ("coronary artery", 1.0)],
    "stroke": [("brain", 1.0)],
    "peripheral arterial disease": [("aorta", 0.8)],   # closest available
    "pulmonary arterial hypertension": [("lung", 1.0), ("aorta", 0.6)],
    "cystic fibrosis": [("lung", 1.0), ("pancreas", 0.6)],
    "polycystic ovary syndrome": [("ovary", 1.0)],
    "osteoporosis": [("bone marrow", 0.6)],
    "gout": [("kidney", 0.5)],
    "hypertension": [("aorta", 0.6), ("kidney", 0.5)],
    "obesity": [("adipose tissue", 1.0)],
    "hypothyroidism": [("thyroid gland", 1.0)],
    "hyperthyroidism": [("thyroid gland", 1.0)],
    "hypercholesterolemia": [("liver", 0.9)],
    "celiac disease": [("small intestine", 1.0), ("duodenum", 1.0)],
    "gastroesophageal reflux disease": [("esophagus", 1.0)],
    "irritable bowel syndrome": [("colon", 1.0), ("small intestine", 0.8)],
    "chronic hepatitis C virus infection": [("liver", 1.0)],
    "hepatitis B virus infection": [("liver", 1.0)],
    "schizophrenia": [("brain", 1.0)],
    "bipolar disorder": [("brain", 1.0)],
    "major depressive disorder": [("brain", 1.0)],
    "generalized anxiety disorder": [("brain", 1.0)],
    "autism spectrum disorder": [("brain", 1.0)],
    "attention deficit hyperactivity disorder": [("brain", 1.0)],
    "sickle cell anemia": [("blood", 1.0), ("bone marrow", 0.7)],
    "duchenne muscular dystrophy": [("skeletal muscle tissue", 1.0), ("heart", 0.8)],
    "spinal muscular atrophy": [("skeletal muscle tissue", 1.0)],
    "marfan syndrome": [("aorta", 1.0), ("heart", 0.7)],
    "ehlers-danlos syndrome": [("skin of body", 0.8)],
    "neurofibromatosis type 1": [("brain", 0.6), ("skin of body", 0.6)],
}


def _parquet_paths(dataset_dir: Path) -> list[str]:
    paths = sorted(dataset_dir.rglob("*.parquet"))
    if not paths:
        raise FileNotFoundError(f"No parquet files in {dataset_dir}")
    return [str(p) for p in paths]


def main() -> None:
    if not DISEASES_DIR.exists():
        raise SystemExit(
            f"{DISEASES_DIR} not found. Run scripts/download_data.sh first.")
    df = pd.read_parquet(_parquet_paths(DISEASES_DIR),
                         columns=["id", "name", "therapeuticAreas"])

    rows = []
    # TA-derived rows
    for _, dis in df.iterrows():
        tas = dis["therapeuticAreas"]
        if tas is None or len(tas) == 0:
            continue
        for ta in tas:
            for tissue, rel in TA_TO_TISSUE.get(ta, []):
                rows.append({"efo_id": dis["id"], "tissue": tissue, "relevance": rel})
    # Disease-name overrides
    df["_name_lc"] = df["name"].astype(str).str.strip().str.lower()
    name_to_id = (df.drop_duplicates(subset=["_name_lc"], keep="first")
                    .set_index("_name_lc")["id"].to_dict())
    for name_lc, tissues in DISEASE_NAME_OVERRIDES.items():
        efo_id = name_to_id.get(name_lc.lower())
        if efo_id is None:
            continue
        for tissue, rel in tissues:
            rows.append({"efo_id": efo_id, "tissue": tissue, "relevance": rel})

    derived = pd.DataFrame(rows)
    # Multiple TAs / overrides for the same disease can produce duplicate
    # (efo_id, tissue) pairs -- keep the MAX relevance per pair.
    derived = (derived.groupby(["efo_id", "tissue"], as_index=False)["relevance"].max())

    if OUT_CSV.exists():
        existing = pd.read_csv(OUT_CSV)
        combined = pd.concat([existing, derived], ignore_index=True)
        combined = combined.drop_duplicates(subset=["efo_id", "tissue"], keep="first")
    else:
        OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
        combined = derived
    combined.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} with {len(combined):,} (efo_id, tissue) rows "
          f"covering {combined['efo_id'].nunique():,} diseases")


if __name__ == "__main__":
    main()
