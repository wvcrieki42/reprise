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
import os
import sqlite3
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable
import pandas as pd


# Env-var fallbacks for secrets. Set these in your shell instead of
# putting tokens in config.full.yaml -- secrets stay out of git, and
# the same config works across dev / prod environments.
NCBI_API_KEY_ENV = "NCBI_API_KEY"
LENS_API_TOKEN_ENV = "LENS_API_TOKEN"


# ----------------------------------------------------------------------
# Defaults -- per-source saturation counts and weights
# ----------------------------------------------------------------------
# At count = saturation, the per-source contribution saturates at 1.0.
DEFAULT_SATURATION = {
    "pubmed": 200,           # abstract-level co-mention (~rare disease ceiling)
    "europepmc": 500,        # full-text catches methods/results mentions
    "clinicaltrials": 10,    # a handful of trials = "this is actively investigated"
    "lens": 50,              # patents are sparser; 50 = heavily protected combo
}
DEFAULT_WEIGHTS = {
    "pubmed": 1.0,
    "europepmc": 1.0,
    "clinicaltrials": 2.0,   # a trial is a stronger signal than a paper
    "lens": 1.5,             # patents are a strong commercial-intent signal
}
DEFAULT_RATE_PER_SEC = {
    "pubmed": 3.0,           # NCBI: 3 req/s without API key, 10 with
    "europepmc": 5.0,
    "clinicaltrials": 5.0,
    "lens": 2.0,             # Lens free tier ~5 req/s; stay polite
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
    lens_api_token: str | None = None      # https://www.lens.org/lens/user/subscriptions
    enable_pubmed: bool = True
    enable_europepmc: bool = True
    enable_clinicaltrials: bool = True
    enable_lens: bool = False              # requires lens_api_token to actually fire
    # Synonym expansion: cap how many alternate terms we OR into each side of the
    # query. PubMed's expression length is finite and 4-each is a sensible cap
    # for recall vs. precision.
    max_synonyms_target: int = 4
    max_synonyms_disease: int = 4
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
        # Env-var fallbacks: explicit config wins, otherwise pick up the
        # shell-exported secret. Empty strings are treated as unset.
        if not self.ncbi_api_key:
            self.ncbi_api_key = os.environ.get(NCBI_API_KEY_ENV) or None
        if not self.lens_api_token:
            self.lens_api_token = os.environ.get(LENS_API_TOKEN_ENV) or None
        # With a token in hand, NCBI lets us go 10 req/s instead of 3.
        if self.ncbi_api_key and self.rate_per_sec.get("pubmed", 0) < 10.0:
            self.rate_per_sec["pubmed"] = 10.0
        # If the operator went to the trouble of setting a Lens token,
        # default-on the patent backend so it actually fires.
        if self.lens_api_token and not self.enable_lens:
            self.enable_lens = True
        self._cache = _CountCache(self.cache_dir / "literature.sqlite",
                                  self.cache_ttl_seconds)
        self._limiters = {s: _RateLimiter(r) for s, r in self.rate_per_sec.items()}
        # session reused across threads (requests.Session is threadsafe for GETs)
        self._session = None

    def credentials_summary(self) -> str:
        """One-line summary of which credentials the client picked up.
        Suitable for logging at the start of a literature pass."""
        bits = []
        bits.append("ncbi_key=" + ("set" if self.ncbi_api_key else "unset"))
        bits.append("lens_token=" + ("set" if self.lens_api_token else "unset"))
        bits.append(f"pubmed_rate={self.rate_per_sec.get('pubmed', 3.0)}/s")
        return ", ".join(bits)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def score_pairs(self, pairs: pd.DataFrame) -> pd.DataFrame:
        """Input columns: target_symbol, efo_id, disease_name.
        Optional: target_synonyms, disease_synonyms (each a semicolon-separated
        string; OR'd into the search query alongside the primary term).
        Output: same keys + pubmed_count, europepmc_count, trial_count,
        patent_count, investigation_prior.
        """
        required = {"target_symbol", "efo_id", "disease_name"}
        missing = required - set(pairs.columns)
        if missing:
            raise KeyError(f"score_pairs: missing columns {missing}")
        keep_cols = ["target_symbol", "efo_id", "disease_name"]
        for opt in ("target_synonyms", "disease_synonyms"):
            if opt in pairs.columns:
                keep_cols.append(opt)
        work = (pairs[keep_cols]
                .dropna(subset=["target_symbol", "efo_id", "disease_name"])
                .drop_duplicates(subset=["target_symbol", "efo_id", "disease_name"])
                .reset_index(drop=True))
        if work.empty:
            return work.assign(pubmed_count=pd.Series(dtype="Int64"),
                               europepmc_count=pd.Series(dtype="Int64"),
                               trial_count=pd.Series(dtype="Int64"),
                               patent_count=pd.Series(dtype="Int64"),
                               investigation_prior=pd.Series(dtype="float"))
        sources_on = [s for s, on in (("pubmed", self.enable_pubmed),
                                      ("europepmc", self.enable_europepmc),
                                      ("clinicaltrials", self.enable_clinicaltrials),
                                      ("lens", self.enable_lens and bool(self.lens_api_token)))
                      if on]
        counts = {s: [None] * len(work) for s in sources_on}
        tasks = []
        for s in sources_on:
            for i, row in work.iterrows():
                t_terms = self._term_list(row.get("target_symbol"),
                                          row.get("target_synonyms", "") if "target_synonyms" in work.columns else "",
                                          self.max_synonyms_target)
                d_terms = self._term_list(row.get("disease_name"),
                                          row.get("disease_synonyms", "") if "disease_synonyms" in work.columns else "",
                                          self.max_synonyms_disease)
                tasks.append((s, i, self._build_query(s, t_terms, d_terms)))
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
        out["patent_count"] = pd.array(counts.get("lens", [pd.NA] * len(work)), dtype="Int64")
        out["investigation_prior"] = [
            self._aggregate({s: counts[s][i] for s in sources_on})
            for i in range(len(work))
        ]
        return out

    # ------------------------------------------------------------------
    # Query construction
    # ------------------------------------------------------------------
    @staticmethod
    def _term_list(primary, synonyms_field, cap: int) -> list[str]:
        """Combine primary term + semicolon-separated synonyms, dedup, cap, sort.

        Sorted so identical (term, synonym) inputs always produce the same query
        string -- which means cache hits are stable across runs.
        """
        primary = "" if primary is None else str(primary).strip()
        terms = [primary] if primary else []
        if synonyms_field:
            for s in str(synonyms_field).split(";"):
                s = s.strip()
                if s and s.lower() not in {t.lower() for t in terms}:
                    terms.append(s)
        # Keep primary first (most distinctive), then alphabetised aliases for stable caching.
        if len(terms) > 1:
            terms = terms[:1] + sorted(terms[1:], key=str.lower)
        return terms[:cap]

    @staticmethod
    def _build_query(source: str, target_terms: list[str], disease_terms: list[str]) -> str:
        if not target_terms or not disease_terms:
            return ""
        if source == "pubmed":
            # tiab = title/abstract -- avoids false hits from author affiliations etc.
            t = " OR ".join(f'"{t}"[tiab]' for t in target_terms)
            d = " OR ".join(f'"{d}"[tiab]' for d in disease_terms)
            return f"({t}) AND ({d})"
        if source == "europepmc":
            t = " OR ".join(f'TITLE_ABS:"{t}"' for t in target_terms)
            d = " OR ".join(f'TITLE_ABS:"{d}"' for d in disease_terms)
            return f"({t}) AND ({d})"
        if source == "clinicaltrials":
            t = " OR ".join(f'"{t}"' for t in target_terms)
            d = " OR ".join(f'"{d}"' for d in disease_terms)
            return f"({t}) AND ({d})"
        if source == "lens":
            # Lens patent search accepts Lucene-style query_string syntax;
            # we search title + abstract + claims for a commercial-intent signal.
            t = " OR ".join(f'"{t}"' for t in target_terms)
            d = " OR ".join(f'"{d}"' for d in disease_terms)
            return f"({t}) AND ({d})"
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
        if source == "lens":
            return self._fetch_lens(query)
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

    def _fetch_lens(self, query: str) -> int:
        """POST to Lens patent search API with a Lucene query_string body.

        Counts patents whose title / abstract / claims co-mention the target
        and disease terms. Requires a personal Lens API token (free tier
        available; set via `lens_api_token`). Returns 0 if unauthorised so
        the run keeps going.
        """
        if not self.lens_api_token:
            return 0
        body = {
            "query": {
                "query_string": {
                    "query": query,
                    "fields": ["title", "abstract", "claims"],
                }
            },
            "size": 0,
        }
        r = self._session.post(
            "https://api.lens.org/patent/search",
            headers={
                "Authorization": f"Bearer {self.lens_api_token}",
                "Content-Type": "application/json",
            },
            json=body, timeout=self.timeout)
        r.raise_for_status()
        return int(r.json().get("total", 0))

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
