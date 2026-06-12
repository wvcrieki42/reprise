"""Populate data/curated/disease_prevalence.csv from a curated prevalence
dictionary + the Open Targets disease name -> EFO ID map.

Why a script:
  Hand-typing EFO IDs for hundreds of diseases is error-prone and breaks
  across OT versions. We keep the source of truth at the disease NAME
  level (lowercase, exact match against OT's canonical name) and look up
  the canonical EFO ID from `data/full/ot_disease_map.csv` (produced by
  scripts/build_full_tables.py).

Behaviour:
  - Existing rows in disease_prevalence.csv are PRESERVED on collision
    (your manual edits win over the script's curated dict).
  - Unmatched names print suggestions from OT (substring match on the
    first word) so you can iterate the dictionary.
  - The OT name match is exact. If a disease is "amyotrophic lateral
    sclerosis" in OT and your key is "ALS", it won't match -- fix the
    key, not the data.

Coverage:
  ~95 diseases of practical interest for US-focused repurposing:
  top SEER cancers, top CDC chronic conditions, common autoimmune /
  neurodegenerative / psychiatric / GI / pulmonary / renal conditions,
  major orphan diseases. Each row cites a source.

Run after:
  bash scripts/download_data.sh
  python scripts/build_full_tables.py
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OT_DISEASE_MAP = ROOT / "data" / "full" / "ot_disease_map.csv"
OUT_CSV = ROOT / "data" / "curated" / "disease_prevalence.csv"

# disease_name (lowercase, exact match against OT) ->
#   (us_patients, us_prevalence_per_100k, source, as_of)
# Prevalence per 100k computed against ~333M US population unless stated
# otherwise. Where official US adult-only figures are used, the rate
# reflects the adult population only; we still report on the full
# denominator for consistency.
CURATED: dict[str, tuple[int, int, str, str]] = {
    # --- Cancers (SEER cancer survivors estimates, 2023) ---
    "breast cancer": (4_100_000, 1235, "SEER cancer survivors estimate", "2023"),
    "prostate cancer": (3_500_000, 1055, "SEER cancer survivors estimate", "2023"),
    "colorectal cancer": (1_500_000, 452, "SEER cancer survivors estimate", "2023"),
    "melanoma": (1_500_000, 452, "SEER cancer survivors estimate", "2023"),
    "lung cancer": (660_000, 199, "SEER cancer survivors estimate", "2023"),
    "thyroid cancer": (1_000_000, 301, "SEER cancer survivors estimate", "2023"),
    "uterine cancer": (870_000, 262, "SEER cancer survivors estimate", "2023"),
    "kidney cancer": (760_000, 229, "SEER cancer survivors estimate", "2023"),
    "urinary bladder cancer": (730_000, 220, "SEER cancer survivors estimate", "2023"),
    "chronic lymphocytic leukemia": (200_000, 60, "SEER / LLS estimates", "2023"),
    "ovarian cancer": (250_000, 75, "SEER cancer survivors estimate", "2023"),
    "pancreatic carcinoma": (100_000, 30, "SEER pancreatic cancer survivors estimate", "2023"),
    "multiple myeloma": (160_000, 48, "SEER / LLS estimates", "2023"),
    "chronic myelogenous leukemia": (70_000, 21, "SEER / LLS estimates", "2023"),
    "acute myeloid leukemia": (90_000, 27, "SEER / LLS estimates", "2023"),
    "acute lymphoblastic leukemia": (60_000, 18, "SEER / LLS estimates", "2023"),
    "gastric cancer": (130_000, 39, "SEER cancer survivors estimate", "2023"),
    "esophageal cancer": (50_000, 15, "SEER cancer survivors estimate", "2023"),
    "hepatocellular carcinoma": (75_000, 23, "SEER liver cancer survivors estimate", "2023"),
    "glioblastoma multiforme": (15_000, 5, "ABTA / NIH brain tumor estimates", "2023"),
    "glioma": (25_000, 8, "ABTA / NIH brain tumor estimates", "2023"),
    "renal cell carcinoma": (760_000, 229, "SEER kidney cancer survivors estimate", "2023"),
    "cervical cancer": (290_000, 87, "SEER cancer survivors estimate", "2023"),
    "head and neck squamous cell carcinoma": (300_000, 90, "SEER oral/pharynx cancer survivors estimate", "2023"),
    "gastrointestinal stromal tumor": (40_000, 12, "GIST patient registry estimates", "2022"),
    # --- Cardiovascular ---
    "hypertension": (116_000_000, 34_980, "CDC NHANES adult prevalence", "2022"),
    "coronary artery disease": (20_000_000, 6_030, "CDC heart disease surveillance", "2022"),
    "heart failure": (6_700_000, 2_020, "AHA heart disease statistics", "2023"),
    "atrial fibrillation": (6_000_000, 1_810, "CDC heart disease surveillance", "2022"),
    "stroke": (7_500_000, 2_260, "CDC stroke surveillance / AHA", "2022"),
    "peripheral arterial disease": (8_500_000, 2_560, "CDC / AHA estimates", "2022"),
    "venous thromboembolism": (900_000, 271, "CDC annual incidence estimate", "2022"),
    # --- Metabolic ---
    "type 1 diabetes mellitus": (2_000_000, 600, "CDC National Diabetes Statistics Report", "2022"),
    "type 2 diabetes mellitus": (37_300_000, 11_230, "CDC National Diabetes Statistics Report", "2022"),
    "diabetes mellitus": (38_400_000, 11_600, "CDC National Diabetes Statistics Report", "2022"),
    "obesity": (108_000_000, 32_500, "CDC NHANES adult obesity prevalence", "2022"),
    "hypercholesterolemia": (28_000_000, 8_430, "CDC heart disease surveillance", "2022"),
    "metabolic syndrome": (66_700_000, 20_090, "CDC NHANES adult prevalence", "2022"),
    # --- Pulmonary ---
    "asthma": (25_000_000, 7_530, "CDC NHIS", "2022"),
    "chronic obstructive pulmonary disease": (16_000_000, 4_820, "CDC BRFSS / NHIS", "2022"),
    "idiopathic pulmonary fibrosis": (100_000, 30, "NIH RareDiseases / PFF estimate", "2020"),
    "pulmonary arterial hypertension": (40_000, 12, "PHA rare-disease estimate", "2022"),
    "cystic fibrosis": (40_000, 12, "CF Foundation patient registry", "2022"),
    "bronchiectasis": (500_000, 151, "American Lung Association estimate", "2022"),
    # --- GI / liver ---
    "inflammatory bowel disease": (1_600_000, 482, "CCFA / CDC estimates", "2020"),
    "crohn's disease": (780_000, 235, "CCFA / NIDDK estimates", "2020"),
    "ulcerative colitis": (907_000, 273, "CCFA / NIDDK estimates", "2020"),
    "irritable bowel syndrome": (30_000_000, 9_040, "Rome Foundation / NIDDK estimates", "2022"),
    "gastroesophageal reflux disease": (60_000_000, 18_080, "NIDDK adult prevalence estimate", "2022"),
    "celiac disease": (3_000_000, 904, "NIDDK / Beyond Celiac estimates", "2022"),
    "non-alcoholic fatty liver disease": (80_000_000, 24_100, "AASLD adult prevalence estimate", "2022"),
    "non-alcoholic steatohepatitis": (5_000_000, 1_510, "AASLD NASH adult prevalence estimate", "2022"),
    "chronic hepatitis c virus infection": (2_400_000, 723, "CDC viral hepatitis surveillance", "2022"),
    "hepatitis b virus infection": (850_000, 256, "CDC viral hepatitis surveillance (chronic carriers)", "2022"),
    # --- Autoimmune / rheumatology ---
    "rheumatoid arthritis": (1_300_000, 392, "CDC / ACR estimates", "2022"),
    "systemic lupus erythematosus": (200_000, 60, "Lupus Foundation estimates", "2022"),
    "psoriasis": (7_500_000, 2_260, "AAD / NPF estimates", "2022"),
    "psoriatic arthritis": (1_500_000, 452, "NPF / ACR estimates", "2022"),
    "ankylosing spondylitis": (400_000, 121, "Spondylitis Association of America estimates", "2022"),
    "sjogren syndrome": (250_000, 75, "Sjogren's Foundation estimates", "2022"),
    "vitiligo": (1_000_000, 301, "AAD estimates", "2022"),
    "alopecia areata": (700_000, 211, "NAAF estimates", "2022"),
    # --- Neurology / psychiatry ---
    "alzheimer disease": (6_700_000, 2_020, "Alzheimer's Association Facts and Figures", "2024"),
    "parkinson disease": (1_000_000, 301, "Parkinson's Foundation estimates", "2023"),
    "huntington disease": (30_000, 9, "Huntington's Disease Society of America", "2022"),
    "amyotrophic lateral sclerosis": (32_000, 10, "ALS Association estimates", "2023"),
    "multiple sclerosis": (1_000_000, 301, "National MS Society estimates", "2022"),
    "lewy body dementia": (1_400_000, 422, "LBDA estimates", "2022"),
    "epilepsy": (3_400_000, 1_025, "CDC adult epilepsy surveillance", "2022"),
    "migraine disorder": (39_000_000, 11_750, "AMF / CDC headache estimates", "2022"),
    "major depressive disorder": (21_000_000, 6_330, "NIMH / SAMHSA estimates", "2022"),
    "schizophrenia": (2_600_000, 783, "NIMH adult prevalence estimate", "2022"),
    "bipolar disorder": (7_000_000, 2_110, "NIMH adult prevalence estimate", "2022"),
    "generalized anxiety disorder": (7_000_000, 2_110, "NIMH adult prevalence estimate", "2022"),
    "autism spectrum disorder": (5_400_000, 1_628, "CDC ADDM Network estimates", "2022"),
    "attention deficit hyperactivity disorder": (10_000_000, 3_010, "CDC NSCH / NIMH estimates", "2022"),
    # --- Renal ---
    "chronic kidney disease": (35_000_000, 10_550, "CDC CKD surveillance system", "2022"),
    "autosomal dominant polycystic kidney disease": (600_000, 181, "PKD Foundation estimates", "2022"),
    "iga glomerulonephritis": (150_000, 45, "NIDDK / NORD estimates (IgA nephropathy)", "2022"),
    # --- Ophthalmology ---
    "age-related macular degeneration": (11_000_000, 3_310, "NEI / AAO estimates", "2022"),
    "glaucoma": (3_000_000, 904, "NEI estimates", "2022"),
    "diabetic retinopathy": (9_600_000, 2_890, "NEI / ADA estimates", "2022"),
    "retinitis pigmentosa": (100_000, 30, "Foundation Fighting Blindness estimates", "2022"),
    # --- Dermatology ---
    "atopic eczema": (32_000_000, 9_640, "AAD / NEA adult+pediatric atopic dermatitis estimates", "2022"),
    "rosacea": (16_000_000, 4_820, "NRS / AAD estimates", "2022"),
    "hidradenitis suppurativa": (350_000, 105, "HS Foundation / NORD estimates", "2022"),
    # --- Infectious ---
    "aids": (1_200_000, 362, "CDC HIV/AIDS surveillance report", "2022"),
    # --- Rare / genetic ---
    "sickle cell anemia": (100_000, 30, "CDC / SCDAA estimates", "2022"),
    "duchenne muscular dystrophy": (15_000, 5, "MDA / NORD estimates", "2022"),
    "spinal muscular atrophy": (25_000, 8, "Cure SMA / MDA estimates", "2022"),
    "hemophilia a": (21_000, 6, "NHF / CDC estimates", "2022"),
    "hemophilia b": (4_000, 1, "NHF / CDC estimates", "2022"),
    "gaucher disease": (6_000, 2, "NGF / NORD estimates", "2022"),
    "fabry disease": (6_000, 2, "NORD / NIH estimates", "2022"),
    "glycogen storage disease ii": (5_000, 2, "AMDA / NORD estimates (Pompe disease)", "2022"),
    "marfan syndrome": (200_000, 60, "Marfan Foundation / NORD estimates", "2022"),
    "ehlers-danlos syndrome": (200_000, 60, "EDS Society estimates", "2022"),
    "neurofibromatosis type 1": (100_000, 30, "CTF / NORD estimates", "2022"),
    # --- Endocrine / other ---
    "hypothyroidism": (16_700_000, 5_030, "ATA adult prevalence estimate", "2022"),
    "hyperthyroidism": (4_700_000, 1_410, "ATA adult prevalence estimate", "2022"),
    "polycystic ovary syndrome": (5_000_000, 1_510, "ENDO Society adult female estimate", "2022"),
    "osteoporosis": (10_000_000, 3_010, "NOF / Bone Health & Osteoporosis Foundation", "2022"),
    "gout": (9_200_000, 2_770, "CDC NHANES / Arthritis Foundation", "2022"),
}


# ----------------------------------------------------------------------
def _id_rank(eid: str) -> tuple[int, int]:
    """EFO_ > MONDO_ > Orphanet_ > DOID_ > HP_ > everything else.
    Within a prefix, lower numeric ID first (older = more foundational)."""
    prefix_order = {"EFO": 0, "MONDO": 1, "Orphanet": 2, "DOID": 3, "HP": 4}
    head, _, tail = eid.partition("_")
    try:
        num = int("".join(c for c in tail if c.isdigit()) or "0")
    except ValueError:
        num = 0
    return (prefix_order.get(head, 99), num)


def main() -> None:
    if not OT_DISEASE_MAP.exists():
        raise SystemExit(
            f"{OT_DISEASE_MAP} not found. Run scripts/download_data.sh and "
            "scripts/build_full_tables.py first."
        )
    ot = pd.read_csv(OT_DISEASE_MAP)
    ot["_name_lc"] = ot["disease_name"].astype(str).str.strip().str.lower()

    rows = []
    unmatched: list[str] = []
    for name_lc, (patients, prev, source, as_of) in CURATED.items():
        hits = ot[ot["_name_lc"] == name_lc]["efo_id"].tolist()
        if not hits:
            unmatched.append(name_lc)
            continue
        efo_id = min(hits, key=_id_rank)
        rows.append({
            "efo_id": efo_id,
            "disease_name": name_lc,
            "us_patients": patients,
            "us_prevalence_per_100k": prev,
            "source": source,
            "as_of": as_of,
        })
    derived = pd.DataFrame(rows)

    # Preserve existing rows (manual edits or earlier-curated entries win).
    if OUT_CSV.exists():
        existing = pd.read_csv(OUT_CSV)
        combined = pd.concat([existing, derived], ignore_index=True)
        combined = combined.drop_duplicates(subset=["efo_id"], keep="first")
    else:
        OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
        combined = derived

    combined.to_csv(OUT_CSV, index=False)

    print(f"matched: {len(rows)} / {len(CURATED)}")
    print(f"wrote {OUT_CSV} with {len(combined)} rows total")
    if unmatched:
        print(f"\nunmatched names ({len(unmatched)}). Suggestions from OT:")
        for u in unmatched:
            first = u.split()[0]
            suggestions = (ot[ot["_name_lc"].str.contains(first, na=False, regex=False)]
                             ["disease_name"].head(5).tolist())
            print(f"  {u!r}")
            for s in suggestions:
                print(f"      -> {s}")


if __name__ == "__main__":
    main()
