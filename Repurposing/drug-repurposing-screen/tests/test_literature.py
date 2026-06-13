"""Unit tests for the LiteraturePriorClient and the two-pass pipeline step.

Uses an in-process fetcher hook so no HTTP is touched.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.config import Config  # noqa: E402
from repurpose.sources.literature import LiteraturePriorClient  # noqa: E402
from repurpose import steps  # noqa: E402


def _client(tmp_path, fetcher, **overrides):
    return LiteraturePriorClient(
        cache_dir=tmp_path,
        fetcher=fetcher,
        cache_ttl_seconds=None,
        max_workers=2,
        **overrides,
    )


def _pairs():
    return pd.DataFrame([
        {"target_symbol": "EGFR",  "efo_id": "EFO_0000305", "disease_name": "breast cancer"},
        {"target_symbol": "OBSCURE", "efo_id": "EFO_9999999", "disease_name": "made up syndrome"},
    ])


def test_aggregation_extremes(tmp_path):
    """count=0 -> prior=0; count at-or-above saturation -> prior=1."""
    saturated = {"EGFR": 10_000, "OBSCURE": 0}

    def fake(source, query):
        return saturated["EGFR"] if "EGFR" in query else saturated["OBSCURE"]

    c = _client(tmp_path, fake)
    out = c.score_pairs(_pairs())
    egfr = out[out.target_symbol == "EGFR"].iloc[0]
    obs = out[out.target_symbol == "OBSCURE"].iloc[0]
    assert egfr["investigation_prior"] == 1.0
    assert obs["investigation_prior"] == 0.0
    # Per-source columns populated; patent_count NaN unless lens is on
    assert int(egfr["pubmed_count"]) == 10_000
    assert int(egfr["europepmc_count"]) == 10_000
    assert int(egfr["trial_count"]) == 10_000
    assert pd.isna(egfr["patent_count"])


def test_disabled_source_omitted(tmp_path):
    calls = []

    def fake(source, query):
        calls.append(source)
        return 100

    c = _client(tmp_path, fake, enable_clinicaltrials=False)
    out = c.score_pairs(_pairs().head(1))
    assert "clinicaltrials" not in calls
    assert pd.isna(out.iloc[0]["trial_count"])
    # Prior averages only over the two enabled sources
    assert out.iloc[0]["investigation_prior"] > 0


def test_cache_reuse(tmp_path):
    calls = []

    def fake(source, query):
        calls.append((source, query))
        return 50

    c1 = _client(tmp_path, fake)
    c1.score_pairs(_pairs().head(1))
    n_first = len(calls)
    assert n_first > 0

    # New client instance, same cache_dir -> no new fetches
    c2 = _client(tmp_path, fake)
    c2.score_pairs(_pairs().head(1))
    assert len(calls) == n_first, "cache should have absorbed all repeat queries"


def test_weighted_aggregation(tmp_path):
    """Trial signal weighted higher than literature signals."""
    def fake(source, query):
        return 10 if source == "clinicaltrials" else 1

    c = _client(tmp_path, fake)
    out = c.score_pairs(_pairs().head(1))
    prior = out.iloc[0]["investigation_prior"]
    # With default weights (pm=1, epmc=1, ct=2) and saturations (200, 500, 10),
    # the trial side dominates and pushes the score well above 0.5.
    assert prior > 0.55


def test_fetch_failure_treated_as_zero(tmp_path):
    def boom(source, query):
        raise RuntimeError("upstream is unhappy")

    c = _client(tmp_path, boom)
    out = c.score_pairs(_pairs().head(1))
    assert out.iloc[0]["investigation_prior"] == 0.0


# ----------------------------------------------------------------------
# Synonym expansion -- target and disease both OR'd into the query
# ----------------------------------------------------------------------
def test_synonym_expansion_in_query(tmp_path):
    """target_synonyms and disease_synonyms get OR'd into the search query."""
    seen_queries = []

    def capture(source, query):
        seen_queries.append((source, query))
        return 1

    c = _client(tmp_path, capture, enable_clinicaltrials=False)
    pairs = pd.DataFrame([{
        "target_symbol": "EGFR",
        "efo_id": "EFO_1",
        "disease_name": "breast cancer",
        "target_synonyms": "Epidermal growth factor receptor;ErbB-1",
        "disease_synonyms": "breast carcinoma;mammary cancer",
    }])
    c.score_pairs(pairs)
    pubmed_q = next(q for s, q in seen_queries if s == "pubmed")
    # All target terms appear OR'd in the target clause
    assert '"EGFR"[tiab]' in pubmed_q
    assert '"Epidermal growth factor receptor"[tiab]' in pubmed_q
    assert '"ErbB-1"[tiab]' in pubmed_q
    # All disease terms appear OR'd in the disease clause
    assert '"breast cancer"[tiab]' in pubmed_q
    assert '"breast carcinoma"[tiab]' in pubmed_q
    assert " OR " in pubmed_q and " AND " in pubmed_q
    # Europe PMC variant uses TITLE_ABS
    epmc_q = next(q for s, q in seen_queries if s == "europepmc")
    assert 'TITLE_ABS:"EGFR"' in epmc_q
    assert 'TITLE_ABS:"breast carcinoma"' in epmc_q


def test_synonym_cap_respected(tmp_path):
    """max_synonyms_target / max_synonyms_disease prevent runaway query bloat."""
    captured = []

    def capture(source, query):
        captured.append(query)
        return 0

    c = _client(tmp_path, capture, enable_europepmc=False, enable_clinicaltrials=False,
                max_synonyms_target=2, max_synonyms_disease=2)
    pairs = pd.DataFrame([{
        "target_symbol": "EGFR",
        "efo_id": "EFO_1",
        "disease_name": "breast cancer",
        "target_synonyms": "alias_a;alias_b;alias_c;alias_d;alias_e",
        "disease_synonyms": "syn_a;syn_b;syn_c;syn_d",
    }])
    c.score_pairs(pairs)
    q = captured[0]
    # Primary + 1 alias on each side = 2 total each; later aliases must be absent
    assert '"alias_a"[tiab]' in q  # alphabetised, kept
    assert '"alias_e"[tiab]' not in q
    assert '"syn_a"[tiab]' in q
    assert '"syn_d"[tiab]' not in q


def test_synonym_ordering_is_stable(tmp_path):
    """Same (primary, synonyms) input -> identical query string -> cache hit."""
    calls = []

    def capture(source, query):
        calls.append((source, query))
        return 5

    c = _client(tmp_path, capture, enable_europepmc=False, enable_clinicaltrials=False)
    base = pd.DataFrame([{
        "target_symbol": "EGFR", "efo_id": "EFO_1", "disease_name": "BC",
        "target_synonyms": "B_alias;A_alias",
    }])
    # Same content, different synonym ordering on second call
    flipped = base.copy()
    flipped.loc[0, "target_synonyms"] = "A_alias;B_alias"
    c.score_pairs(base)
    c.score_pairs(flipped)
    # Either second call hits the cache (no new fetcher invocations) or the
    # query strings are identical -- both prove sort-then-cap is deterministic.
    pubmed_calls = [q for s, q in calls if s == "pubmed"]
    assert len(set(pubmed_calls)) == 1, f"unstable query strings: {pubmed_calls}"


# ----------------------------------------------------------------------
# Lens patent backend
# ----------------------------------------------------------------------
def test_lens_off_by_default(tmp_path):
    """Without lens_api_token, lens isn't queried even if enable_lens=True."""
    seen_sources = set()

    def capture(source, query):
        seen_sources.add(source)
        return 1

    # enable_lens=True but no token -> still skipped
    c = _client(tmp_path, capture, enable_lens=True)
    c.score_pairs(_pairs().head(1))
    assert "lens" not in seen_sources


def test_lens_fires_when_token_set(tmp_path):
    """With token + enable_lens, lens is queried and patent_count is populated."""
    def fake(source, query):
        if source == "lens":
            return 25
        return 5

    c = _client(tmp_path, fake, enable_lens=True, lens_api_token="fake-token-for-test")
    out = c.score_pairs(_pairs().head(1))
    assert int(out.iloc[0]["patent_count"]) == 25
    assert out.iloc[0]["investigation_prior"] > 0


# ----------------------------------------------------------------------
# Credentials -- env-var fallback + auto-enable-on-token behaviour
# ----------------------------------------------------------------------
from unittest.mock import patch as _patch


def test_env_var_fallback_for_ncbi_key(tmp_path):
    """No config-level key + NCBI_API_KEY in env -> picked up automatically."""
    with _patch.dict("os.environ", {"NCBI_API_KEY": "env-ncbi-key-123",
                                    "LENS_API_TOKEN": ""}, clear=False):
        c = LiteraturePriorClient(cache_dir=tmp_path)
        assert c.ncbi_api_key == "env-ncbi-key-123"
        # And throughput should bump to 10 req/s
        assert c.rate_per_sec["pubmed"] == 10.0


def test_explicit_config_overrides_env(tmp_path):
    """Explicit constructor arg wins over env var."""
    with _patch.dict("os.environ", {"NCBI_API_KEY": "env-key"}, clear=False):
        c = LiteraturePriorClient(cache_dir=tmp_path, ncbi_api_key="explicit-key")
        assert c.ncbi_api_key == "explicit-key"


def test_env_var_fallback_for_lens_token(tmp_path):
    """LENS_API_TOKEN in env -> token picked up AND enable_lens flipped on."""
    with _patch.dict("os.environ", {"LENS_API_TOKEN": "env-lens-token",
                                    "NCBI_API_KEY": ""}, clear=False):
        c = LiteraturePriorClient(cache_dir=tmp_path)
        assert c.lens_api_token == "env-lens-token"
        assert c.enable_lens is True


def test_credentials_summary(tmp_path):
    c = LiteraturePriorClient(cache_dir=tmp_path,
                              ncbi_api_key="k", lens_api_token="t")
    s = c.credentials_summary()
    assert "ncbi_key=set" in s
    assert "lens_token=set" in s
    assert "pubmed_rate=10.0" in s


def test_no_credentials_keeps_defaults(tmp_path):
    """No env, no config -> nothing set, throughput stays at default."""
    with _patch.dict("os.environ", {"NCBI_API_KEY": "", "LENS_API_TOKEN": ""},
                     clear=False):
        c = LiteraturePriorClient(cache_dir=tmp_path)
        assert c.ncbi_api_key is None
        assert c.lens_api_token is None
        assert c.enable_lens is False
        assert c.rate_per_sec["pubmed"] == 3.0


# ----------------------------------------------------------------------
# Two-pass step integration -- exercises steps.add_literature_pass
# ----------------------------------------------------------------------
def _ranked():
    """Synthetic ranked hypotheses, already sorted by opportunity desc."""
    return pd.DataFrame([
        # rank 1-3: top-K territory; rank 4-5: outside the top-K we'll set
        {"rank": 1, "drug_id": "D1", "efo_id": "E1", "disease_name": "well known disease",
         "lead_target": "TGT_HOT", "opportunity": 0.90},
        {"rank": 2, "drug_id": "D2", "efo_id": "E2", "disease_name": "novel disease",
         "lead_target": "TGT_COLD", "opportunity": 0.80},
        {"rank": 3, "drug_id": "D3", "efo_id": "E3", "disease_name": "medium disease",
         "lead_target": "TGT_MID", "opportunity": 0.70},
        {"rank": 4, "drug_id": "D4", "efo_id": "E4", "disease_name": "out of band",
         "lead_target": "TGT_X", "opportunity": 0.30},
        {"rank": 5, "drug_id": "D5", "efo_id": "E5", "disease_name": "way out",
         "lead_target": "TGT_Y", "opportunity": 0.20},
    ])


def _cfg(top_n=3, w_inv=0.8):
    return Config(
        raw={"literature": {"top_n": top_n}, "scoring": {"w_investigation": w_inv}},
        root=Path("."),
    )


def test_two_pass_only_queries_top_k(tmp_path):
    queried_targets = set()

    def fake(source, query):
        if "TGT_HOT" in query:
            queried_targets.add("TGT_HOT")
            return 500
        if "TGT_COLD" in query:
            queried_targets.add("TGT_COLD")
            return 0
        if "TGT_MID" in query:
            queried_targets.add("TGT_MID")
            return 20
        # Anything outside top-K must NOT be queried
        queried_targets.add("LEAKED")
        return 999

    client = _client(tmp_path, fake)
    out = steps.add_literature_pass(_ranked(), client, _cfg(top_n=3))
    assert queried_targets == {"TGT_HOT", "TGT_COLD", "TGT_MID"}, \
        f"queried beyond top-K: {queried_targets}"
    # Rows outside top-K have no literature data
    out_of_band = out[out.drug_id.isin(["D4", "D5"])]
    assert out_of_band["investigation_prior"].isna().all()
    assert out_of_band["opportunity"].tolist() == [0.30, 0.20]


def test_two_pass_reorders_by_damping(tmp_path):
    """The well-investigated row gets pushed below the genuinely novel one."""
    def fake(source, query):
        if "TGT_HOT" in query:
            return 10_000   # saturated -> investigation_prior near 1
        if "TGT_COLD" in query:
            return 0
        return 5            # weak

    client = _client(tmp_path, fake)
    out = steps.add_literature_pass(_ranked(), client, _cfg(top_n=3, w_inv=0.9))
    # TGT_HOT/D1 started rank 1 with opp 0.90 but should fall behind TGT_COLD/D2
    top_drug = out.iloc[0]["drug_id"]
    assert top_drug == "D2", f"expected D2 (novel) on top, got {top_drug}"
    # And its rank is 1 after re-ranking
    assert int(out.iloc[0]["rank"]) == 1
    # D1 opportunity got damped substantially
    d1 = out[out.drug_id == "D1"].iloc[0]
    assert d1["opportunity"] < 0.20, f"expected heavy damping, got {d1['opportunity']}"
    assert d1["investigation_prior"] > 0.95


def test_two_pass_zero_weight_is_inspection_only(tmp_path):
    """w_investigation = 0 -> attach columns but don't change opportunity / order."""
    def fake(source, query):
        return 1000 if "TGT_HOT" in query else 0

    client = _client(tmp_path, fake)
    original = _ranked()
    out = steps.add_literature_pass(original, client, _cfg(top_n=3, w_inv=0.0))
    # Order unchanged
    assert out["drug_id"].tolist() == original["drug_id"].tolist()
    # Opportunity unchanged
    assert (out["opportunity"].values == original["opportunity"].values).all()
    # But literature columns are populated for the top-K
    top_k = out[out.drug_id.isin(["D1", "D2", "D3"])]
    assert top_k["investigation_prior"].notna().all()


# ----------------------------------------------------------------------
# Extended Lens-only coverage (lens_top_n > top_n)
# ----------------------------------------------------------------------
def test_only_sources_restricts_to_named_sources(tmp_path):
    """score_pairs(only_sources={'lens'}) hits only the Lens fetcher and
    populates only patent_count + investigation_prior."""
    seen = []

    def fake(source, query):
        seen.append(source)
        return 7

    c = _client(tmp_path, fake, enable_lens=True, lens_api_token="t")
    out = c.score_pairs(_pairs().head(1), only_sources={"lens"})
    assert set(seen) == {"lens"}, f"expected lens-only, got {set(seen)}"
    row = out.iloc[0]
    assert int(row["patent_count"]) == 7
    assert pd.isna(row["pubmed_count"])
    assert pd.isna(row["europepmc_count"])
    assert pd.isna(row["trial_count"])
    # investigation_prior is computed from the Lens count alone
    assert 0 < row["investigation_prior"] <= 1


def _cfg_with_lens_extension(top_n=3, lens_top_n=5, w_inv=0.8):
    return Config(
        raw={"literature": {"top_n": top_n, "lens_top_n": lens_top_n},
             "scoring": {"w_investigation": w_inv}},
        root=Path("."),
    )


def test_lens_top_n_extends_patent_coverage(tmp_path):
    """With lens_top_n > top_n, Lens fires for rows beyond top_n but the
    cheap sources (PubMed / EPMC / NCT) do NOT."""
    calls = []

    def fake(source, query):
        calls.append((source, query))
        # any non-zero count is fine for the assertion
        return 50

    client = _client(tmp_path, fake, enable_lens=True,
                     lens_api_token="fake-token-for-test")
    out = steps.add_literature_pass(_ranked(), client,
                                    _cfg_with_lens_extension(top_n=3,
                                                             lens_top_n=5,
                                                             w_inv=0.0))
    sources_by_target = {}
    for source, query in calls:
        for t in ("TGT_HOT", "TGT_COLD", "TGT_MID", "TGT_X", "TGT_Y"):
            if t in query:
                sources_by_target.setdefault(t, set()).add(source)
    # Top-3 targets get all four sources
    for t in ("TGT_HOT", "TGT_COLD", "TGT_MID"):
        assert {"pubmed", "europepmc", "clinicaltrials", "lens"} <= sources_by_target[t], \
            f"{t} only saw {sources_by_target.get(t)}"
    # Extended targets get LENS ONLY -- the cheap sources stay capped at top_n
    for t in ("TGT_X", "TGT_Y"):
        assert sources_by_target[t] == {"lens"}, \
            f"{t} should be lens-only past top_n, got {sources_by_target.get(t)}"
    # And patent_count is populated for the extended rows
    ext = out[out.drug_id.isin(["D4", "D5"])]
    assert ext["patent_count"].notna().all()
    # investigation_prior also populated for the extended rows (Lens signal only)
    assert ext["investigation_prior"].notna().all()


def test_lens_top_n_default_equals_top_n(tmp_path):
    """Omitting lens_top_n -> behaves identically to the legacy single-pass."""
    seen_targets = set()

    def fake(source, query):
        for t in ("TGT_HOT", "TGT_COLD", "TGT_MID", "TGT_X", "TGT_Y"):
            if t in query:
                seen_targets.add(t)
        return 10

    client = _client(tmp_path, fake, enable_lens=True, lens_api_token="t")
    steps.add_literature_pass(_ranked(), client, _cfg(top_n=3))
    assert seen_targets == {"TGT_HOT", "TGT_COLD", "TGT_MID"}, \
        f"Lens leaked past top_n when lens_top_n was unset: {seen_targets}"


if __name__ == "__main__":
    import tempfile
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            with tempfile.TemporaryDirectory() as d:
                fn(Path(d))
                print(f"OK {name}")
