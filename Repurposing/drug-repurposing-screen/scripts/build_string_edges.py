"""Build symbol-keyed STRING edge list from the downloaded protein.links file.

STRING ships protein-protein interaction scores keyed on Ensembl protein
IDs (9606.ENSP00000xxxxxx). The repurposing pipeline expects HGNC
symbols in `source_symbol,target_symbol,score` form. This script does
the conversion:

  1. Read 9606.protein.info.v12.0.txt.gz for ENSP -> preferred_name
     (HGNC symbol) mapping. Downloads it on first run if not present.
  2. Stream 9606.protein.links.v12.0.txt.gz, filter score >= MIN_SCORE,
     translate both endpoints to symbols, emit one row per kept edge.

Output: data/full/string_edges.csv with columns source_symbol,
target_symbol, score (score normalised to [0, 1] from STRING's 0-1000).

Run after scripts/download_data.sh. Idempotent -- rerun to refresh.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FULL = ROOT / "data" / "full"
LINKS_PATH = FULL / "9606.protein.links.v12.0.txt.gz"
INFO_PATH = FULL / "9606.protein.info.v12.0.txt.gz"
OUT_PATH = FULL / "string_edges.csv"
INFO_URL = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"


def _ensure_info_file() -> None:
    if INFO_PATH.exists():
        return
    print(f"downloading {INFO_URL} -> {INFO_PATH}")
    INFO_PATH.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(INFO_URL, INFO_PATH)


def _load_ensp_to_symbol() -> dict[str, str]:
    """Parse the STRING info file: ENSP id -> preferred_name (HGNC symbol)."""
    mapping: dict[str, str] = {}
    with gzip.open(INFO_PATH, "rt") as fh:
        header = fh.readline()
        # columns: string_protein_id preferred_name protein_size annotation
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            ensp, sym = parts[0], parts[1]
            # Strip the 9606. species prefix so it matches the links file
            mapping[ensp] = sym
    return mapping


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-score", type=int, default=400,
                    help="STRING combined_score (0-1000) cutoff; 400 = medium, "
                         "700 = high (default 400 = pipeline-side filter does the rest)")
    ap.add_argument("--output", default=str(OUT_PATH))
    args = ap.parse_args()

    if not LINKS_PATH.exists():
        raise SystemExit(
            f"{LINKS_PATH} not found. Run scripts/download_data.sh first.")
    _ensure_info_file()

    print(f"loading ENSP -> symbol from {INFO_PATH.name}")
    ensp_to_sym = _load_ensp_to_symbol()
    print(f"  mapped {len(ensp_to_sym):,} ENSP ids -> symbol")

    print(f"streaming {LINKS_PATH.name} (min_score={args.min_score})")
    n_in, n_out, n_unmapped = 0, 0, 0
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(LINKS_PATH, "rt") as fh, open(out_path, "w", newline="") as out:
        header = fh.readline()  # "protein1 protein2 combined_score"
        writer = csv.writer(out)
        writer.writerow(["source_symbol", "target_symbol", "score"])
        for line in fh:
            n_in += 1
            parts = line.rstrip("\n").split(" ")
            if len(parts) != 3:
                continue
            p1, p2, s = parts
            score = int(s)
            if score < args.min_score:
                continue
            sym1 = ensp_to_sym.get(p1)
            sym2 = ensp_to_sym.get(p2)
            if not sym1 or not sym2:
                n_unmapped += 1
                continue
            if sym1 == sym2:
                continue
            writer.writerow([sym1, sym2, round(score / 1000.0, 3)])
            n_out += 1
    print(f"  read {n_in:,} edges, wrote {n_out:,}; {n_unmapped:,} unmapped "
          f"({n_unmapped / max(n_in, 1) * 100:.1f}%)")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
