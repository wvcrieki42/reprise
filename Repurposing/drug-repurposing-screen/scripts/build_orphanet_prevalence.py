"""Build data/curated/orphanet_prevalence.csv from Orphanet's bulk prevalence XML.

For each disorder in en_product9_prev.xml, pick the best prevalence entry
and convert the prevalence-class range to an approximate US patient count.
The output CSV has the same schema as data/curated/disease_prevalence.csv
so the MarketSizeClient can layer it underneath the manually curated rows.

Best-entry priority (per disorder):
  1. Geographic: United States > Worldwide > Europe > other
  2. Type: Point prevalence > Annual incidence > Birth prevalence > Lifetime > Cases/families
  3. Class is one of the structured ranges (excludes None / Unknown / Not yet documented)

PrevalenceClass -> midpoint per-100k population (assuming ~333M US):

    >1 / 1000          -> 100 per 100k    (~333k US patients; floor)
    1-5 / 10 000       ->  30 per 100k    (~100k)
    6-9 / 10 000       ->  75 per 100k    (~250k)
    1-9 / 100 000      ->   5 per 100k    (~16.7k)
    1-9 / 1 000 000    ->   0.5 per 100k  (~1.7k)
    <1 / 1 000 000     ->   0.05 per 100k (~167 -- ultra-rare orphan)

Run after `bash scripts/download_data.sh` -- expects
data/full/orphanet/en_product9_prev.xml on disk.
"""
from __future__ import annotations
import csv
import glob
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XML_PATH = ROOT / "data" / "full" / "orphanet" / "en_product9_prev.xml"
OT_DISEASES_DIR = ROOT / "data" / "full" / "diseases"
OUT_CSV = ROOT / "data" / "curated" / "orphanet_prevalence.csv"

US_POPULATION = 333_000_000

# Geographic preference: higher score wins
GEO_PRIORITY = {
    "United States of America": 4,
    "United States": 4,
    "Worldwide": 3,
    "Europe": 2,
    # any others -> 1
}

# Type preference: higher score wins
TYPE_PRIORITY = {
    "Point prevalence": 5,
    "Annual incidence": 4,
    "Prevalence at birth": 3,
    "Lifetime Prevalence": 2,
    "Cases/families": 1,
}

# PrevalenceClass -> per-100k midpoint
CLASS_TO_PER100K = {
    ">1 / 1000": 100.0,             # >100 per 100k -- use floor
    "1-5 / 10 000": 30.0,           # midpoint of [10, 50]
    "6-9 / 10 000": 75.0,           # midpoint of [60, 90]
    "1-9 / 100 000": 5.0,           # midpoint of [1, 9]
    "1-9 / 1 000 000": 0.5,         # midpoint of [0.1, 0.9]
    "<1 / 1 000 000": 0.05,         # ultra-rare
}


def _text(node):
    return node.text if node is not None else None


def _entry_score(geo: str | None, type_: str | None) -> tuple[int, int]:
    return (GEO_PRIORITY.get(geo, 1), TYPE_PRIORITY.get(type_, 0))


def _build_orpha_to_ot_id() -> dict[str, str]:
    """Map Orpha code -> OT primary disease id (usually MONDO_xxx) via OT's dbXRefs.

    The pipeline keys diseases by OT's primary id; Orphanet's bulk file
    uses Orpha codes (e.g. 166024). The Orphanet -> OT bridge lives in
    OT's diseases parquet under the `dbXRefs` field.
    """
    if not OT_DISEASES_DIR.exists():
        return {}
    import pandas as pd
    files = sorted(glob.glob(str(OT_DISEASES_DIR / "*.parquet")))
    if not files:
        return {}
    df = pd.read_parquet(files, columns=["id", "dbXRefs"])
    mapping: dict[str, str] = {}
    for ot_id, xrefs in zip(df["id"], df["dbXRefs"]):
        if xrefs is None:
            continue
        for x in xrefs:
            s = str(x)
            lower = s.lower()
            if lower.startswith("orphanet:") or lower.startswith("orphanet_"):
                code = s.split(":", 1)[-1].split("_", 1)[-1].strip()
                # First write wins -- OT often lists the same Orpha against
                # both an EFO and a MONDO id; we want MONDO/EFO over DOID etc.
                # but a simple first-wins is good enough at this scale.
                if code and code not in mapping:
                    mapping[code] = ot_id
    return mapping


def main() -> None:
    if not XML_PATH.exists():
        raise SystemExit(f"{XML_PATH} not found. Run bash scripts/download_data.sh first.")
    root = ET.parse(XML_PATH).getroot()
    orpha_to_ot = _build_orpha_to_ot_id()
    print(f"loaded {len(orpha_to_ot):,} Orphanet -> OT id mappings from "
          f"{OT_DISEASES_DIR}")

    rows = []
    n_mapped, n_unmapped = 0, 0
    skipped = 0
    for disorder in root.iter("Disorder"):
        orpha = _text(disorder.find("OrphaCode"))
        name = _text(disorder.find("Name"))
        if not orpha or not name:
            continue
        best_entry = None
        best_score = (-1, -1)
        for prev in disorder.findall("PrevalenceList/Prevalence"):
            geo = _text(prev.find("PrevalenceGeographic/Name"))
            type_ = _text(prev.find("PrevalenceType/Name"))
            cls = _text(prev.find("PrevalenceClass/Name"))
            if cls not in CLASS_TO_PER100K:
                continue
            score = _entry_score(geo, type_)
            if score > best_score:
                best_score = score
                best_entry = (geo, type_, cls)
        if best_entry is None:
            skipped += 1
            continue
        geo, type_, cls = best_entry
        per_100k = CLASS_TO_PER100K[cls]
        us_patients = int(round(per_100k / 100_000 * US_POPULATION))
        # Prefer the OT primary id (MONDO_xxx) when we have it -- pipeline
        # diseases are keyed that way. Fall back to Orphanet_xxx so we
        # don't drop rows the bridge doesn't cover (~1k rows).
        ot_id = orpha_to_ot.get(orpha)
        if ot_id is not None:
            efo_id = ot_id
            n_mapped += 1
        else:
            efo_id = f"Orphanet_{orpha}"
            n_unmapped += 1
        rows.append({
            "efo_id": efo_id,
            "disease_name": name,
            "us_patients": us_patients,
            "us_prevalence_per_100k": per_100k,
            "source": f"Orphanet {type_} ({cls}, {geo})",
            "as_of": "2025",
        })

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "efo_id", "disease_name", "us_patients",
            "us_prevalence_per_100k", "source", "as_of"])
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {OUT_CSV} with {len(rows):,} rows "
          f"(skipped {skipped:,} disorders without a usable prevalence class)")
    print(f"  keyed by OT primary id (MONDO/EFO): {n_mapped:,}")
    print(f"  keyed by Orphanet_xxx fallback:     {n_unmapped:,}")


if __name__ == "__main__":
    main()
