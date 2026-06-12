"""Command-line entry point.

Usage:
    python -m repurpose.cli --config config.yaml [--top 50] [--quiet]
"""
from __future__ import annotations
import argparse
from .pipeline import run_from_file


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Mechanism-driven drug-repurposing screen")
    ap.add_argument("--config", default="config.yaml", help="path to config.yaml")
    ap.add_argument("--top", type=int, default=25, help="rows to preview to stdout")
    ap.add_argument("--quiet", action="store_true", help="suppress progress logging")
    args = ap.parse_args(argv)

    ranked = run_from_file(args.config, verbose=not args.quiet)
    if len(ranked):
        import pandas as pd
        pd.set_option("display.max_colwidth", 40)
        pd.set_option("display.width", 200)
        preview = ranked.head(args.top)[
            ["rank", "drug_name", "disease_name", "mechanistic_support",
             "novelty", "opportunity", "evidence_targets"]]
        print("\nTop hypotheses:\n")
        print(preview.to_string(index=False))
    else:
        print("No hypotheses passed the thresholds.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
