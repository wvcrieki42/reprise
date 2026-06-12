"""STRING functional-interaction network expansion.

Two backends:
  * file : read a pre-downloaded edge list (offline / reproducible / demo)
  * api  : query the live STRING REST API, with on-disk caching

Edge-list schema (file backend): source_symbol, target_symbol, score   (score 0..1)
"""
from __future__ import annotations
from pathlib import Path
import json
import pandas as pd


class StringClient:
    def __init__(self, source: str = "file", *, edges_path: Path | None = None,
                 min_confidence: float = 0.7, max_partners: int = 10,
                 species: int = 9606, cache_dir: Path | None = None):
        self.source = source
        self.min_confidence = float(min_confidence)
        self.max_partners = int(max_partners)
        self.species = species
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self._edges: pd.DataFrame | None = None
        if source == "file":
            if not edges_path or not Path(edges_path).exists():
                # network expansion gracefully degrades to "no partners"
                self._edges = pd.DataFrame(columns=["source_symbol", "target_symbol", "score"])
            else:
                e = pd.read_csv(edges_path)
                e["score"] = pd.to_numeric(e["score"], errors="coerce").fillna(0.0)
                self._edges = e
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    def neighbors(self, symbols: list[str]) -> pd.DataFrame:
        """Return DataFrame[target_symbol, partner_symbol, score] for the inputs."""
        symbols = sorted({s for s in symbols if s})
        if not symbols:
            return pd.DataFrame(columns=["target_symbol", "partner_symbol", "score"])
        if self.source == "file":
            return self._neighbors_file(symbols)
        return self._neighbors_api(symbols)

    def _neighbors_file(self, symbols: list[str]) -> pd.DataFrame:
        e = self._edges
        m = e[(e["source_symbol"].isin(symbols)) & (e["score"] >= self.min_confidence)]
        m = m.rename(columns={"source_symbol": "target_symbol", "target_symbol": "partner_symbol"})
        m = m[m["partner_symbol"] != m["target_symbol"]]
        return (m.sort_values("score", ascending=False)
                 .groupby("target_symbol", group_keys=False)
                 .head(self.max_partners)[["target_symbol", "partner_symbol", "score"]])

    def _neighbors_api(self, symbols: list[str]) -> pd.DataFrame:
        import requests  # local import so demo mode needs no network deps
        key = None
        if self.cache_dir:
            key = self.cache_dir / (json.dumps(symbols)[:200].replace("/", "_") + ".json")
            if key.exists():
                return pd.read_json(key)
        ids = "%0d".join(symbols)
        url = (f"https://string-db.org/api/json/interaction_partners?identifiers={ids}"
               f"&species={self.species}&required_score={int(self.min_confidence*1000)}"
               f"&limit={self.max_partners}&caller_identity=repurpose_pipeline")
        rows = requests.get(url, timeout=60).json()
        out = pd.DataFrame([{
            "target_symbol": r["preferredName_A"],
            "partner_symbol": r["preferredName_B"],
            "score": r["score"],
        } for r in rows])
        if out.empty:
            out = pd.DataFrame(columns=["target_symbol", "partner_symbol", "score"])
        if key is not None:
            out.to_json(key)
        return out
