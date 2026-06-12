"""KOL finder: one US + one EU key-opinion leader per (target, disease) pair.

For each pair, query PubMed for recent papers, aggregate authors across them,
then pick the highest-publishing author per region (US, EU) with affiliation
preference for:

  US: Boston / Cambridge MA / Bay Area (SF, Stanford, Berkeley, Palo Alto)
  EU: Belgium > Netherlands > France / Germany / Luxembourg > UK > other Europe

Best author affiliation -> institution + email. h-index from Semantic
Scholar's free Author API.

Per-pair cost: 1 esearch + 1 efetch (PubMed) + ~2 Semantic Scholar lookups
(one per region). With an NCBI API key (10 req/s) and Semantic Scholar's
generous public quota, top_n=100 finishes in ~20-30s; top_n=1000 in a
few minutes. SQLite cache makes re-runs near-instant.
"""
from __future__ import annotations
import math
import os
import re
import sqlite3
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable
import pandas as pd


NCBI_API_KEY_ENV = "NCBI_API_KEY"
SEMANTIC_SCHOLAR_API_KEY_ENV = "SEMANTIC_SCHOLAR_API_KEY"


# ----------------------------------------------------------------------
# Region classification + affiliation scoring
# ----------------------------------------------------------------------
_EU_COUNTRIES = {
    "BELGIUM", "BELGIQUE", "BELGIË",
    "NETHERLANDS", "HOLLAND",
    "FRANCE", "GERMANY", "DEUTSCHLAND", "LUXEMBOURG",
    "UNITED KINGDOM", "UK", "ENGLAND", "SCOTLAND", "WALES",
    "IRELAND", "SPAIN", "ITALY", "PORTUGAL",
    "AUSTRIA", "SWITZERLAND", "DENMARK", "SWEDEN", "NORWAY", "FINLAND",
    "POLAND", "CZECH REPUBLIC", "CZECHIA", "HUNGARY", "GREECE",
    "ROMANIA", "BULGARIA", "CROATIA", "SLOVAKIA", "SLOVENIA",
}
_US_NAMES = {"USA", "U.S.A.", "U.S.A", "U.S.", "UNITED STATES", "UNITED STATES OF AMERICA"}

_US_CITY_BONUS: dict[str, int] = {  # higher = better
    # Boston cluster
    "BOSTON": 100, "CAMBRIDGE": 95,            # MA only -- guarded below
    "HARVARD": 100, "MIT": 100, "MGH": 90, "MASSACHUSETTS GENERAL": 100,
    "DANA-FARBER": 95, "DANA FARBER": 95, "BRIGHAM": 90, "TUFTS": 80,
    "BROAD INSTITUTE": 100, "WHITEHEAD": 90, "BETH ISRAEL DEACONESS": 85,
    "CHILDREN'S HOSPITAL BOSTON": 95, "BOSTON CHILDREN": 95,
    # Bay Area cluster
    "SAN FRANCISCO": 100, "STANFORD": 100, "BERKELEY": 100, "PALO ALTO": 95,
    "MOUNTAIN VIEW": 80, "UCSF": 100, "GLADSTONE": 90, "BUCK INSTITUTE": 80,
}
_EU_CITY_BONUS: dict[str, int] = {
    # Belgium
    "BELGIUM": 100, "BELGIQUE": 100, "BELGIË": 100,
    "GHENT": 100, "GENT": 100, "UGENT": 100, "KU LEUVEN": 100, "LEUVEN": 100,
    "BRUSSELS": 100, "BRUXELLES": 100, "ULB": 100, "VUB": 100, "UCL": 100,
    "ANTWERP": 100, "ANTWERPEN": 100, "LIÈGE": 100, "LIEGE": 100,
    # Netherlands (closest neighbour)
    "NETHERLANDS": 80, "AMSTERDAM": 80, "ROTTERDAM": 80, "UTRECHT": 80,
    "MAASTRICHT": 80, "LEIDEN": 80, "NIJMEGEN": 80, "GRONINGEN": 80,
    # France / Germany / Luxembourg
    "PARIS": 70, "STRASBOURG": 70, "LILLE": 70, "LYON": 70,
    "GERMANY": 70, "DEUTSCHLAND": 70, "BERLIN": 70, "MUNICH": 70,
    "HEIDELBERG": 70, "AACHEN": 70, "COLOGNE": 70, "KÖLN": 70,
    "LUXEMBOURG": 70,
    # UK
    "LONDON": 60, "OXFORD": 60, "CAMBRIDGE UK": 60, "EDINBURGH": 60,
}

# Anything in PubMed that looks like a CONTACT email
_EMAIL_RE = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")


def _classify_region(affiliation: str) -> str:
    """Return 'US' / 'EU' / '' for a PubMed affiliation string."""
    if not affiliation:
        return ""
    up = affiliation.upper()
    # Guard for "Cambridge MA" vs "Cambridge UK"
    if ", MA" in up or "MASSACHUSETTS" in up:
        return "US"
    # Country-based first
    for country in _US_NAMES:
        if country in up:
            return "US"
    for country in _EU_COUNTRIES:
        if country in up:
            return "EU"
    return ""


def _affiliation_score(affiliation: str, region: str) -> int:
    """Higher = more preferred (Boston/Bay for US, Belgium-cluster for EU)."""
    if not affiliation:
        return 0
    up = affiliation.upper()
    table = _US_CITY_BONUS if region == "US" else _EU_CITY_BONUS
    score = 0
    for needle, weight in table.items():
        if needle in up:
            score = max(score, weight)
    return score


def _extract_email(affiliation: str) -> str:
    if not affiliation:
        return ""
    m = _EMAIL_RE.search(affiliation)
    return m.group(0) if m else ""


def _shorten_institution(affiliation: str) -> str:
    """First clause of the affiliation, stripped of email."""
    if not affiliation:
        return ""
    cleaned = _EMAIL_RE.sub("", affiliation).strip().strip(".")
    # PubMed convention: comma-separated, first piece is the department / first
    # is institution. Heuristic: take the first chunk that looks like an org.
    parts = [p.strip() for p in cleaned.split(",") if p.strip()]
    # Drop bare-department-style leading parts ("Department of X") -- keep the
    # first part that mentions a university / institute / hospital / etc.
    org_hints = ("UNIVERSITY", "INSTITUTE", "HOSPITAL", "CENTER", "CENTRE",
                 "SCHOOL", "COLLEGE", "FACULTY", "LABORATORY", "RESEARCH")
    for p in parts:
        if any(h in p.upper() for h in org_hints):
            return p[:160]
    return parts[0][:160] if parts else cleaned[:160]


def _match_scholar_authors(wanted_name: str, items: list[dict]) -> list[tuple[int, int]]:
    """Score Semantic Scholar author results against a PubMed name.

    PubMed: 'Eric Topol' / 'E Topol' / 'Pieter Janssens'
    SS:     'E. Topol'  / 'Eric Topol' / 'P. Janssens' / 'E.J. Topol'

    Match rule: same last name (exact), and first-name first-initial agrees.
    Score 100 when full first name matches, 80 when only initial matches.
    """
    parts = wanted_name.strip().split()
    if not parts:
        return []
    last = parts[-1].lower()
    first = parts[0].lower() if parts else ""
    first_initial = first[0] if first else ""

    scored: list[tuple[int, int]] = []
    for item in items:
        cand_name = (item.get("name") or "").strip()
        cand_parts = cand_name.split()
        if not cand_parts:
            continue
        cand_last = cand_parts[-1].lower()
        if cand_last != last:
            continue
        # Examples to handle for the first-name token:
        #   'Eric' vs 'Eric'         -> full match
        #   'Eric' vs 'E.' / 'E'     -> initial match
        #   'E' vs 'Eric'            -> initial match
        cand_first = cand_parts[0].lower().rstrip(".")
        score = 0
        if cand_first == first and len(first) > 1:
            score = 100
        elif first_initial and cand_first and cand_first[0] == first_initial:
            score = 80
        if score > 0:
            scored.append((score, item.get("hIndex") or 0))
    # Both "Boris Keren" and "B. Keren" can be the same person -- and when the
    # name is ambiguous (4-5 candidates with the same surname + initial),
    # the most-cited author IS overwhelmingly the more likely KOL match.
    # Sort by h-index DESC primarily, name-quality only as a tiebreak.
    scored.sort(key=lambda s: (s[1], s[0]), reverse=True)
    return scored


# ----------------------------------------------------------------------
# Cache
# ----------------------------------------------------------------------
class _Cache:
    def __init__(self, path: Path, ttl_seconds: float | None):
        self.path = Path(path)
        self.ttl = ttl_seconds
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        with self._connect() as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS kol_results (
                  query_key  TEXT PRIMARY KEY,
                  payload    TEXT NOT NULL,
                  fetched_at REAL NOT NULL
                )
            """)

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.path, timeout=30.0)

    def get(self, key: str) -> str | None:
        with self._lock, self._connect() as conn:
            row = conn.execute(
                "SELECT payload, fetched_at FROM kol_results WHERE query_key = ?",
                (key,)).fetchone()
        if row is None:
            return None
        payload, fetched_at = row
        if self.ttl is not None and (time.time() - fetched_at) > self.ttl:
            return None
        return payload

    def put(self, key: str, payload: str) -> None:
        with self._lock, self._connect() as conn:
            conn.execute(
                "INSERT OR REPLACE INTO kol_results (query_key, payload, fetched_at) "
                "VALUES (?, ?, ?)",
                (key, payload, time.time()))


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
# Client
# ----------------------------------------------------------------------
@dataclass
class KOLClient:
    """One US + one EU KOL per (target, disease) query, with h-index when possible."""
    cache_dir: Path
    ncbi_api_key: str | None = None
    semantic_scholar_api_key: str | None = None
    max_pmids_per_pair: int = 50
    max_workers: int = 4
    timeout: float = 30.0
    cache_ttl_seconds: float | None = 90 * 24 * 3600
    # Test hooks
    pubmed_fetcher: Callable[[str], str] | None = None
    semantic_scholar_fetcher: Callable[[str], dict] | None = None

    def __post_init__(self) -> None:
        if not self.ncbi_api_key:
            self.ncbi_api_key = os.environ.get(NCBI_API_KEY_ENV) or None
        if not self.semantic_scholar_api_key:
            self.semantic_scholar_api_key = os.environ.get(SEMANTIC_SCHOLAR_API_KEY_ENV) or None
        self.cache_dir = Path(self.cache_dir)
        self._cache = _Cache(self.cache_dir / "kol.sqlite", self.cache_ttl_seconds)
        # NCBI: 3 req/s without key, 10 with
        ncbi_rate = 10.0 if self.ncbi_api_key else 3.0
        self._ncbi_limiter = _RateLimiter(ncbi_rate)
        self._scholar_limiter = _RateLimiter(1.0)  # Semantic Scholar polite default
        self._session = None

    # ------------------------------------------------------------------
    def find_kols(self, pairs: pd.DataFrame) -> pd.DataFrame:
        """Input cols: target_symbol, efo_id, disease_name, [disease_synonyms].
        Output: same keys + us_kol_name/institution/email/h_index/n_pubs +
        the same five fields prefixed eu_."""
        required = {"target_symbol", "efo_id", "disease_name"}
        missing = required - set(pairs.columns)
        if missing:
            raise KeyError(f"find_kols: missing columns {missing}")
        keep = ["target_symbol", "efo_id", "disease_name"]
        if "disease_synonyms" in pairs.columns:
            keep.append("disease_synonyms")
        work = (pairs[keep]
                .dropna(subset=["target_symbol", "efo_id", "disease_name"])
                .drop_duplicates(subset=["target_symbol", "efo_id", "disease_name"])
                .reset_index(drop=True))
        if work.empty:
            return work.assign(**{f"{r}_kol_{c}": pd.Series(dtype=object)
                                  for r in ("us", "eu")
                                  for c in ("name", "institution", "email")})

        results = [None] * len(work)
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futures = {pool.submit(self._kol_for_pair,
                                   row.target_symbol, row.disease_name,
                                   row.get("disease_synonyms", "") if "disease_synonyms" in work.columns else ""): i
                       for i, row in work.iterrows()}
            for fut, i in futures.items():
                try:
                    results[i] = fut.result()
                except Exception:
                    results[i] = {"us": {}, "eu": {}}

        out = work.copy()
        for region in ("us", "eu"):
            out[f"{region}_kol_name"] = [r[region].get("name", "") for r in results]
            out[f"{region}_kol_institution"] = [r[region].get("institution", "") for r in results]
            out[f"{region}_kol_email"] = [r[region].get("email", "") for r in results]
            out[f"{region}_kol_h_index"] = pd.array(
                [r[region].get("h_index") for r in results], dtype="Int64")
            out[f"{region}_kol_n_pubs"] = pd.array(
                [r[region].get("n_pubs", 0) for r in results], dtype="Int64")
        return out

    # ------------------------------------------------------------------
    def _kol_for_pair(self, target: str, disease: str,
                      disease_synonyms: str = "") -> dict:
        key = f"{target}||{disease}||{disease_synonyms or ''}"
        cached = self._cache.get(key)
        if cached is not None:
            import json
            return json.loads(cached)
        # Build query
        query = self._build_query(target, disease, disease_synonyms)
        pmids = self._pubmed_esearch(query)
        if not pmids:
            payload = {"us": {}, "eu": {}}
        else:
            articles = self._pubmed_efetch(pmids[:self.max_pmids_per_pair])
            payload = self._select_kols(articles)
        import json
        self._cache.put(key, json.dumps(payload))
        return payload

    @staticmethod
    def _build_query(target: str, disease: str, synonyms: str = "") -> str:
        d_terms = [disease]
        if synonyms:
            for s in synonyms.split(";"):
                s = s.strip()
                if s and s.lower() != disease.lower():
                    d_terms.append(s)
        d_clause = " OR ".join(f'"{t}"[tiab]' for t in d_terms[:4])
        return f'("{target}"[tiab]) AND ({d_clause})'

    # ------------------------------------------------------------------
    # PubMed
    # ------------------------------------------------------------------
    def _pubmed_esearch(self, query: str) -> list[str]:
        params = {"db": "pubmed", "term": query, "retmode": "json",
                  "retmax": self.max_pmids_per_pair, "sort": "relevance"}
        if self.ncbi_api_key:
            params["api_key"] = self.ncbi_api_key
        self._ncbi_limiter.wait()
        if self.pubmed_fetcher is not None:
            data = self.pubmed_fetcher("esearch:" + query)
            if isinstance(data, str):
                import json
                data = json.loads(data)
        else:
            import requests
            if self._session is None:
                self._session = requests.Session()
            r = self._session.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params=params, timeout=self.timeout)
            r.raise_for_status()
            data = r.json()
        return list(data.get("esearchresult", {}).get("idlist", []))

    def _pubmed_efetch(self, pmids: list[str]) -> list[dict]:
        if not pmids:
            return []
        self._ncbi_limiter.wait()
        if self.pubmed_fetcher is not None:
            xml_str = self.pubmed_fetcher("efetch:" + ",".join(pmids))
        else:
            import requests
            if self._session is None:
                self._session = requests.Session()
            params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
            if self.ncbi_api_key:
                params["api_key"] = self.ncbi_api_key
            r = self._session.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params=params, timeout=self.timeout)
            r.raise_for_status()
            xml_str = r.text
        return self._parse_efetch_xml(xml_str)

    @staticmethod
    def _parse_efetch_xml(xml_str: str) -> list[dict]:
        articles: list[dict] = []
        if not xml_str:
            return articles
        try:
            root = ET.fromstring(xml_str)
        except ET.ParseError:
            return articles
        for art in root.iter("PubmedArticle"):
            authors: list[dict] = []
            for auth in art.iter("Author"):
                last = auth.findtext("LastName") or ""
                fore = auth.findtext("ForeName") or auth.findtext("Initials") or ""
                if not last:
                    continue
                affs = [a.text for a in auth.iter("Affiliation") if a.text]
                authors.append({
                    "name": (f"{fore} {last}").strip(),
                    "affiliation": " | ".join(affs),
                })
            year = art.findtext(".//PubDate/Year") or ""
            articles.append({"authors": authors, "year": year})
        return articles

    # ------------------------------------------------------------------
    # KOL selection
    # ------------------------------------------------------------------
    def _select_kols(self, articles: list[dict]) -> dict:
        """Aggregate authors across articles, pick the highest-scoring per region."""
        # author_name -> {n_pubs, best_affiliation, best_affiliation_score, region}
        agg: dict[str, dict] = {}
        for art in articles:
            for a in art["authors"]:
                name = a["name"]
                if not name:
                    continue
                affiliation = a["affiliation"] or ""
                region = _classify_region(affiliation)
                if not region:
                    continue
                score = _affiliation_score(affiliation, region)
                rec = agg.setdefault(name, {
                    "n_pubs": 0, "best_aff": "", "best_score": -1, "region": ""
                })
                rec["n_pubs"] += 1
                if score > rec["best_score"]:
                    rec["best_score"] = score
                    rec["best_aff"] = affiliation
                    rec["region"] = region

        # Pick the best author per region.
        # Primary sort: affiliation score (region preference). Tiebreak: n_pubs.
        # Tiebreak again: name (stable).
        chosen: dict[str, dict] = {"us": {}, "eu": {}}
        for region in ("us", "eu"):
            candidates = [(name, rec) for name, rec in agg.items()
                          if rec["region"] == region.upper()]
            if not candidates:
                continue
            name, rec = max(candidates,
                            key=lambda x: (x[1]["best_score"], x[1]["n_pubs"], -len(x[0])))
            picked = {
                "name": name,
                "institution": _shorten_institution(rec["best_aff"]),
                "email": _extract_email(rec["best_aff"]),
                "n_pubs": rec["n_pubs"],
            }
            h = self._semantic_scholar_h_index(name)
            if h is not None:
                picked["h_index"] = h
            chosen[region] = picked
        return chosen

    # ------------------------------------------------------------------
    # Semantic Scholar h-index
    # ------------------------------------------------------------------
    def _semantic_scholar_h_index(self, name: str) -> int | None:
        if not name:
            return None
        self._scholar_limiter.wait()
        if self.semantic_scholar_fetcher is not None:
            data = self.semantic_scholar_fetcher(name)
        else:
            import requests
            if self._session is None:
                self._session = requests.Session()
            headers = {}
            if self.semantic_scholar_api_key:
                headers["x-api-key"] = self.semantic_scholar_api_key
            try:
                r = self._session.get(
                    "https://api.semanticscholar.org/graph/v1/author/search",
                    params={"query": name, "limit": 3, "fields": "name,hIndex"},
                    headers=headers, timeout=self.timeout)
                r.raise_for_status()
                data = r.json()
            except Exception:
                return None
        if not isinstance(data, dict):
            return None
        items = data.get("data", [])
        if not items:
            return None
        scored = _match_scholar_authors(name, items)
        if not scored or scored[0][0] < 60:
            return None
        h = scored[0][1]
        return int(h) if h else None
