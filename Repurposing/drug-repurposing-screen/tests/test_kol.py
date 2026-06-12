"""Unit tests for the KOL (key-opinion-leader) finder.

Mocks the PubMed efetch + Semantic Scholar API so no HTTP is touched.
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.sources.kol import (  # noqa: E402
    KOLClient, _classify_region, _affiliation_score, _extract_email,
    _shorten_institution,
)


# ---------------------------------------------------------------- helpers
def test_classify_region_us_canonical():
    assert _classify_region("Massachusetts General Hospital, Boston, MA, USA") == "US"
    assert _classify_region("Stanford University, Palo Alto, CA, United States") == "US"


def test_classify_region_eu_belgium():
    assert _classify_region("Ghent University, Ghent, Belgium") == "EU"
    assert _classify_region("KU Leuven, Belgique") == "EU"


def test_classify_region_cambridge_disambiguation():
    """'Cambridge, MA' must be US; 'Cambridge, UK' must be EU."""
    assert _classify_region("MIT, Cambridge, MA 02139, USA") == "US"
    assert _classify_region("University of Cambridge, Cambridge, UK") == "EU"


def test_classify_region_unknown():
    assert _classify_region("Some private clinic, Tokyo, Japan") == ""
    assert _classify_region("") == ""


def test_affiliation_score_prefers_boston_for_us():
    boston = _affiliation_score("Harvard Medical School, Boston, MA, USA", "US")
    other_us = _affiliation_score("University of Iowa, Iowa City, IA, USA", "US")
    assert boston > other_us


def test_affiliation_score_prefers_belgium_for_eu():
    belgium = _affiliation_score("KU Leuven, Belgium", "EU")
    italy = _affiliation_score("University of Rome, Rome, Italy", "EU")
    assert belgium > italy


def test_email_extraction():
    aff = "Dana-Farber Cancer Institute, Boston, MA, USA. jane.smith@dfci.harvard.edu"
    assert _extract_email(aff) == "jane.smith@dfci.harvard.edu"
    assert _extract_email("No email here") == ""


def test_shorten_institution_finds_organisation_chunk():
    aff = "Department of Pediatrics, Massachusetts General Hospital, Boston, MA, USA"
    assert "Massachusetts General Hospital" in _shorten_institution(aff)


# ---------------------------------------------------------------- end-to-end with stub
def _stub_pubmed_factory(articles_xml: str):
    """Returns a fetcher that handles esearch and efetch URLs."""
    def fetch(url_key: str) -> str:
        if url_key.startswith("esearch:"):
            return json.dumps({"esearchresult": {"idlist": ["1", "2", "3"]}})
        if url_key.startswith("efetch:"):
            return articles_xml
        return ""
    return fetch


def _stub_scholar_factory(h_index_by_name: dict):
    def fetch(name: str) -> dict:
        h = h_index_by_name.get(name)
        if h is None:
            return {"data": []}
        return {"data": [{"name": name, "hIndex": h}]}
    return fetch


def _articles_xml() -> str:
    """3 articles: 2 by a Boston author + 1 by a Leuven author."""
    return """<?xml version="1.0"?>
<PubmedArticleSet>
  <PubmedArticle><MedlineCitation><Article>
    <AuthorList>
      <Author>
        <LastName>Smith</LastName>
        <ForeName>Jane</ForeName>
        <AffiliationInfo><Affiliation>Dana-Farber Cancer Institute, Boston, MA 02115, USA. jane.smith@dfci.harvard.edu</Affiliation></AffiliationInfo>
      </Author>
    </AuthorList>
    <PubDate><Year>2024</Year></PubDate>
  </Article></MedlineCitation></PubmedArticle>
  <PubmedArticle><MedlineCitation><Article>
    <AuthorList>
      <Author>
        <LastName>Smith</LastName>
        <ForeName>Jane</ForeName>
        <AffiliationInfo><Affiliation>Dana-Farber Cancer Institute, Boston, MA, USA</Affiliation></AffiliationInfo>
      </Author>
    </AuthorList>
    <PubDate><Year>2023</Year></PubDate>
  </Article></MedlineCitation></PubmedArticle>
  <PubmedArticle><MedlineCitation><Article>
    <AuthorList>
      <Author>
        <LastName>Janssens</LastName>
        <ForeName>Pieter</ForeName>
        <AffiliationInfo><Affiliation>KU Leuven, Belgium. pieter.janssens@kuleuven.be</Affiliation></AffiliationInfo>
      </Author>
    </AuthorList>
    <PubDate><Year>2024</Year></PubDate>
  </Article></MedlineCitation></PubmedArticle>
</PubmedArticleSet>"""


def test_end_to_end_finds_one_us_and_one_eu(tmp_path):
    client = KOLClient(
        cache_dir=tmp_path,
        pubmed_fetcher=_stub_pubmed_factory(_articles_xml()),
        semantic_scholar_fetcher=_stub_scholar_factory(
            {"Jane Smith": 42, "Pieter Janssens": 27}),
    )
    pairs = pd.DataFrame([{
        "target_symbol": "EGFR",
        "efo_id": "EFO_X",
        "disease_name": "breast cancer",
    }])
    out = client.find_kols(pairs)
    row = out.iloc[0]
    assert row["us_kol_name"] == "Jane Smith"
    assert "Dana-Farber" in row["us_kol_institution"]
    assert row["us_kol_email"] == "jane.smith@dfci.harvard.edu"
    assert int(row["us_kol_h_index"]) == 42
    assert int(row["us_kol_n_pubs"]) == 2          # 2 articles by Smith
    assert row["eu_kol_name"] == "Pieter Janssens"
    assert "KU Leuven" in row["eu_kol_institution"]
    assert row["eu_kol_email"] == "pieter.janssens@kuleuven.be"
    assert int(row["eu_kol_h_index"]) == 27


def test_cache_reuse(tmp_path):
    calls = []
    def stub_pubmed(key):
        calls.append(key)
        if key.startswith("esearch"):
            return json.dumps({"esearchresult": {"idlist": ["1"]}})
        return _articles_xml()
    client = KOLClient(
        cache_dir=tmp_path, pubmed_fetcher=stub_pubmed,
        semantic_scholar_fetcher=_stub_scholar_factory({"Jane Smith": 42}),
    )
    pairs = pd.DataFrame([{"target_symbol": "EGFR", "efo_id": "EFO_X",
                           "disease_name": "breast cancer"}])
    client.find_kols(pairs)
    n_first = len(calls)
    # Second instance reading from same cache_dir -> no new PubMed calls
    client2 = KOLClient(
        cache_dir=tmp_path, pubmed_fetcher=stub_pubmed,
        semantic_scholar_fetcher=_stub_scholar_factory({"Jane Smith": 42}),
    )
    client2.find_kols(pairs)
    assert len(calls) == n_first


def test_no_kol_when_no_us_or_eu_authors(tmp_path):
    xml = """<?xml version="1.0"?>
<PubmedArticleSet><PubmedArticle><MedlineCitation><Article>
  <AuthorList><Author>
    <LastName>Tanaka</LastName><ForeName>Hiroshi</ForeName>
    <AffiliationInfo><Affiliation>University of Tokyo, Tokyo, Japan</Affiliation></AffiliationInfo>
  </Author></AuthorList>
  <PubDate><Year>2024</Year></PubDate>
</Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"""
    client = KOLClient(
        cache_dir=tmp_path,
        pubmed_fetcher=_stub_pubmed_factory(xml),
        semantic_scholar_fetcher=lambda n: {"data": []},
    )
    out = client.find_kols(pd.DataFrame([{
        "target_symbol": "EGFR", "efo_id": "EFO_X", "disease_name": "breast cancer"
    }]))
    row = out.iloc[0]
    assert row["us_kol_name"] == ""
    assert row["eu_kol_name"] == ""


if __name__ == "__main__":
    import tempfile
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            sig = fn.__code__.co_varnames[:fn.__code__.co_argcount]
            with tempfile.TemporaryDirectory() as d:
                if "tmp_path" in sig:
                    fn(Path(d))
                else:
                    fn()
                print(f"OK {name}")
