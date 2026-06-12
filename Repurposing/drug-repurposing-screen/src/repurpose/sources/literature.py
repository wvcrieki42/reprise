"""Literature / patent / trial prior for (target, disease) pairs.

Answers the question the current `novelty` score completely misses:
*"has anyone actually investigated this target for this disease?"*
A pair can be ontologically novel (no approved indication, no near-neighbour
in EFO) yet have 500 PubMed hits and three Phase 2 trials. We want to know.

The client queries three free, no-auth APIs and aggregates the hit counts:

  * PubMed E-utilities  -- abstract-level co-mention count
  * Europe PMC           -- includes OA full-text co-mentions
  * ClinicalTrials.gov   -- trial-level investigation (much stronger signal)

Each source is normalised to [0, 1] via a log-saturating function and the
weighted mean is the per-pair `investigation_prior` (1 = heavily studied,
0 = no prior mention). Lens / patent sources can be plugged in later via
the same `_fetch_*` -> `(int)` contract.

Cache: a SQLite database under `cache_dir/literature.sqlite`. Each
(source, query) tuple is fetched at most once; subsequent runs reuse the
count. Counts on these sources change slowly enough that a 90-day TTL is
generous (configurable).

Designed to be called once per *unique* (target_symbol, efo_id) pair --
the same target reappears across many drugs, so this cache is the main
reason the pipeline stays tractable.
"""
from __future__ import annotations
import math
import sqlite3
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable
import pandas as pd


# ----------------------------------------------------------------------
# Defaults -- per-source saturation counts and weights
# ----------------------------------------------------------------------
# At count = saturation, the per-source contribution saturates at 1.0.
DEFAULT_SATURATION = {
    "pubmed": 200,           # abstract-level co-mention (~rare disease ceiling)
    "europepmc": 500,        # full-text catches methods/results mentions
    "clinicaltrials": 10,    # a handful of trials = "this is actively investigated"
}
DEFAULT_WEIGHTS = {
    "pubmed": 1.0,
    "europepmc": 1.0,
    "clinicaltrials": 2.0,   # a trial is a stronger signal than a paper
}
DEFAULT_RATE_PER_SEC = {
    "pubmed": 3.0,           # NCBI: 3 req/s without API key, 10 with
    "europepmc": 5.0,
    "clinicaltrials": 5.0,
}


# ----------------------------------------------------------------------
# Rate limiter -- simple token bucket per source
# ----------------------------------------------------------------------
class _RateLimiter:
    def __init__(self, rate_per_sec: float):
        self._interval = 1.0 / max(rate_per_sec, 0.1)
        self._next_ok = 0.0
        self._lock = threading.Lock()

    def wait(self) -> None:
        with self._lock:
            now = time.monotonic()
            sleep_for = self._next_ok - now
            if sleep_for > 0:
                time.sleep(sleep_for)
                now += sleep_for
            self._next_ok = now + self._interval


# ----------------------------------------------------------------------
# Cache -- SQLite, (source, query) -> count
# ----------------------------------------------------------------------
class _CountCache:
    def __init__(self, path: Path, ttl_seconds: float | None):
        self.path = Path(path)
        self.ttl = ttl_seconds
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        with self._connect() as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS counts (
                  source     TEXT NOT NULL,
                  query      TEXT NOT NULL,
                  count      INTEGER NOT NULL,
                  fetched_at REAL NOT NULL,
                  PRIMARY KEY (source, query)
                )
            """)

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.path, timeout=30.0)

    def get(self, source: str, query: str) -> int | None:
        with self._lock, self._connect() as conn:
            row = conn.execute(
                "SELECT count, fetched_at FROM counts WHERE source = ? AND query = ?",
                (source, query),
            ).fetchone()
        if row is None:
            return None
        count, fetched_at = row
        if self.ttl is not None and (time.time() - fetched_at) > self.ttl:
            return None
        return int(count)

    def put(self, source: str, query: str, count: int) -> None:
        with self._lock, self._connect() as conn:
            conn.execute(
                "INSERT OR REPLACE INTO counts (source, query, count, fetched_at) VALUES (?, ?, ?, ?)",
                (source, query, int(count), time.time()),
            )


# ----------------------------------------------------------------------
# Client
# ----------------------------------------------------------------------
@dataclass
class LiteraturePriorClient:
    """Score (target, disease) pairs by literature / patent / trial evidence.

    Drop the client into the pipeline, hand it a DataFrame of unique
    (target_symbol, efo_id, disease_name) rows, and get back the same
    rows plus per-source counts and an aggregated `investigation_prior`.
    """
    cache_dir: Path
    ncbi_api_key: str | None = None
    enable_pubmed: bool = True
    enable_europepmc: bool = True
    enable_clinicaltrials: bool = True
    saturation: dict[str, int] = field(default_factory=lambda: dict(DEFAULT_SATURATION))
    weights: dict[str, float] = field(default_factory=lambda: dict(DEFAULT_WEIGHTS))
    rate_per_sec: dict[str, float] = field(default_factory=lambda: dict(DEFAULT_RATE_PER_SEC))
    max_workers: int = 4
    timeout: float = 20.0
    cache_ttl_seconds: float | None = 90 * 24 * 3600  # 90 days
    # Test hook: swap in a stub `(source, query) -> int` to avoid HTTP in tests.
    fetcher: Callable[[str, str], int] | None = None

    def __post_init__(self) -> None:
        self.cache_dir = Path(self.cache_dir)
        self._cache = _CountCache(self.cache_dir / "literature.sqlite",
                                  self.cache_ttl_seconds)
        self._limiters = {s: _RateLimiter(r) for s, r in self.rate_per_sec.items()}
        # session reused across threads (requests.Session is threadsafe for GETs)
        self._session = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def score_pairs(self, pairs: pd.DataFrame) -> pd.DataFrame:
        """Input columns: target_symbol, efo_id, disease_name.
        Output: same keys + pubmed_count, europepmc_count, trial_count,
        investigation_prior (all per-source counts default to NaN if disabled).
        """
        required = {"target_symbol", "efo_id", "disease_name"}
        missing = required - set(pairs.columns)
        if missing:
            raise KeyError(f"score_pairs: missing columns {missing}")
        work = (pairs[["target_symbol", "efo_id", "disease_name"]]
                .dropna()
                .drop_duplicates()
                .reset_index(drop=True))
        if work.empty:
            return work.assign(pubmed_count=pd.Series(dtype="Int64"),
                               europepmc_count=pd.Series(dtype="Int64"),
                               trial_count=pd.Series(dtype="Int64"),
                               investigation_prior=pd.Series(dtype="float"))
        sources_on = [s for s, on in (("pubmed", self.enable_pubmed),
                                      ("europepmc", self.enable_europepmc),
                                      ("clinicaltrials", self.enable_clinicaltrials)) if on]
        counts = {s: [None] * len(work) for s in sources_on}
        tasks = [(s, i, self._build_query(s, row.target_symbol, row.disease_name))
                 for s in sources_on
                 for i, row in work.iterrows()]
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futures = {pool.submit(self._fetch_cached, s, q): (s, i)
                       for s, i, q in tasks}
            for fut in futures:
                s, i = futures[fut]
                try:
                    counts[s][i] = int(fut.result())
                except Exception:
                    counts[s][i] = 0   # treat fetch failure as "no evidence", don't crash a run
        out = work.copy()
        out["pubmed_count"] = pd.array(counts.get("pubmed", [pd.NA] * len(work)), dtype="Int64")
        out["europepmc_count"] = pd.array(counts.get("europepmc", [pd.NA] * len(work)), dtype="Int64")
        out["trial_count"] = pd.array(counts.get("clinicaltrials", [pd.NA] * len(work)), dtype="Int64")
        out["investigation_prior"] = [
            self._aggregate({s: counts[s][i] for s in sources_on})
            for i in range(len(work))
        ]
        return out

    # ------------------------------------------------------------------
    # Query construction
    # ------------------------------------------------------------------
    @staticmethod
    def _build_query(source: str, target: str, disease: str) -> str:
        t = str(target).strip()
        d = str(disease).strip()
        if source == "pubmed":
            # tiab = title/abstract -- avoids false hits from author affiliations etc.
            return f'("{t}"[tiab]) AND ("{d}"[tiab])'
        if source == "europepmc":
            return f'(TITLE_ABS:"{t}" AND TITLE_ABS:"{d}")'
        if source == "clinicaltrials":
            # ClinicalTrials.gov v2 free-text query; AND-combines automatically
            return f'"{t}" AND "{d}"'
        raise ValueError(f"unknown source: {source}")

    # ------------------------------------------------------------------
    # Cached fetch dispatcher
    # ------------------------------------------------------------------
    def _fetch_cached(self, source: str, query: str) -> int:
        cached = self._cache.get(source, query)
        if cached is not None:
            return cached
        self._limiters[source].wait()
        if self.fetcher is not None:
            count = int(self.fetcher(source, query))
        else:
            count = int(self._fetch_live(source, query))
        self._cache.put(source, query, count)
        return count

    def _fetch_live(self, source: str, query: str) -> int:
        # Local import so test runs and the offline demo don't require `requests`.
        import requests
        if self._session is None:
            self._session = requests.Session()
        if source == "pubmed":
            return self._fetch_pubmed(query)
        if source == "europepmc":
            return self._fetch_europepmc(query)
        if source == "clinicaltrials":
            return self._fetch_clinicaltrials(query)
        raise ValueError(f"unknown source: {source}")

    def _fetch_pubmed(self, query: str) -> int:
        params = {"db": "pubmed", "term": query, "rettype": "count", "retmode": "json"}
        if self.ncbi_api_key:
            params["api_key"] = self.ncbi_api_key
        r = self._session.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params=params, timeout=self.timeout)
        r.raise_for_status()
        return int(r.json().get("esearchresult", {}).get("count", 0))

    def _fetch_europepmc(self, query: str) -> int:
        r = self._session.get(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            params={"query": query, "format": "json", "pageSize": 1},
            timeout=self.timeout)
        r.raise_for_status()
        return int(r.json().get("hitCount", 0))

    def _fetch_clinicaltrials(self, query: str) -> int:
        r = self._session.get(
            "https://clinicaltrials.gov/api/v2/studies",
            params={"query.term": query, "pageSize": 1, "countTotal": "true"},
            timeout=self.timeout)
        r.raise_for_status()
        return int(r.json().get("totalCount", 0))

    # ------------------------------------------------------------------
    # Aggregation
    # ------------------------------------------------------------------
    def _aggregate(self, counts: dict[str, int | None]) -> float:
        num, denom = 0.0, 0.0
        for source, count in counts.items():
            if count is None:
                continue
            sat = max(int(self.saturation.get(source, 100)), 1)
            norm = min(math.log1p(max(count, 0)) / math.log1p(sat), 1.0)
            w = float(self.weights.get(source, 1.0))
            num += w * norm
            denom += w
        if denom == 0.0:
            return 0.0
        return round(num / denom, 4)
