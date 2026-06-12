"""US market-size lookup for predicted diseases.

Inspection-only by design: attaches per-disease US patient counts and
prevalence-per-100k to ranked hypotheses without biasing the opportunity
score. Lets you sort/filter hypotheses by addressable market alongside
mechanistic_support and novelty.

V1 backend: a curated CSV at `cfg.path("disease_prevalence")` keyed by
efo_id. The CSV is the single source of truth -- you populate it from
SEER / CDC / Orphanet / GBD / hand-curated estimates and the pipeline
just joins on efo_id. Diseases without an entry get NaN (preserved
through to the output) -- we don't fabricate numbers we don't have.

The orphan-drug threshold (US patients < 200,000 = FDA Orphan Drug Act
definition) is a separate inspection-only flag. Use it to surface
small-but-strategically-valuable opportunities (orphan exclusivity is a
real market lever) without forcing them into the ranking score.

Future backends sketched in `backends`:
  * wikidata    : SPARQL P2854 (prevalence) / P1603 (cases) lookup
  * gbd         : IHME GHDx API (requires registration token)
  * orphanet    : rare-disease prevalence (bulk XML)
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable
import pandas as pd


CURATED_REQUIRED = {"efo_id", "us_patients", "us_prevalence_per_100k", "source", "as_of"}


@dataclass
class MarketSizeClient:
    """Looks up per-disease US market data. V1: curated CSV only."""
    curated_csv: Path | None = None
    backends: list[str] = field(default_factory=lambda: ["curated_csv"])

    def __post_init__(self) -> None:
        self._csv: pd.DataFrame = pd.DataFrame(columns=list(CURATED_REQUIRED))
        if "curated_csv" in self.backends and self.curated_csv:
            p = Path(self.curated_csv)
            if p.exists():
                df = pd.read_csv(p)
                missing = CURATED_REQUIRED - set(df.columns)
                if missing:
                    raise KeyError(
                        f"market CSV missing required columns {missing}: {p}")
                df["us_patients"] = pd.to_numeric(df["us_patients"], errors="coerce")
                df["us_prevalence_per_100k"] = pd.to_numeric(
                    df["us_prevalence_per_100k"], errors="coerce")
                df["source"] = df["source"].fillna("").astype(str)
                df["as_of"] = df["as_of"].fillna("").astype(str)
                self._csv = df.drop_duplicates(subset=["efo_id"], keep="first")

    # ------------------------------------------------------------------
    def lookup(self, efo_ids: Iterable[str]) -> pd.DataFrame:
        """Return per-EFO US market data for the requested ids.

        Columns: efo_id, us_patients (Int64, NA on miss), us_prevalence_per_100k
        (float, NaN on miss), market_source (str), as_of (str).
        Diseases without a CSV entry come back with NaN/empty values --
        callers can decide how to display that ("unknown" vs. "0").
        """
        ids = pd.Series(list(efo_ids), name="efo_id").dropna().drop_duplicates().reset_index(drop=True)
        base = pd.DataFrame({"efo_id": ids})
        if self._csv.empty:
            base["us_patients"] = pd.array([pd.NA] * len(base), dtype="Int64")
            base["us_prevalence_per_100k"] = pd.array(
                [float("nan")] * len(base), dtype="float")
            base["market_source"] = ""
            base["as_of"] = ""
            return base
        csv = self._csv.rename(columns={"source": "market_source"})
        merged = base.merge(
            csv[["efo_id", "us_patients", "us_prevalence_per_100k", "market_source", "as_of"]],
            on="efo_id", how="left")
        merged["us_patients"] = merged["us_patients"].astype("Int64")
        merged["market_source"] = merged["market_source"].fillna("")
        merged["as_of"] = merged["as_of"].fillna("")
        return merged
