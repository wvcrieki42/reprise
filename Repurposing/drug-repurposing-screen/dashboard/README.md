# REPRISE dashboard

**Live: https://reprise.streamlit.app**

Public Streamlit companion to the REPRISE drug-repurposing screen. Three tabs:

1. **Browse hypotheses** -- the 176,272-row ranked output, sidebar filters,
   opportunity-vs-mech scatter, click any point or row for the full per-hit
   brief (mechanism evidence + IP runway + literature & patent prior +
   clinical opportunity + chemical structure from ChEMBL + PubMed search
   link + downloadable PDF).
2. **History & how REPRISE adds value** -- timeline of the 54 historic
   repurposing wins used to validate the screen, plus a disease-disease
   network where edges are the drugs that bridge two indications.
3. **FAQ** -- one expander per computed column with the actual formula,
   honest caveats, and a one-paragraph paper TL;DR.

## Run locally

```bash
pip install -r dashboard/requirements.txt
streamlit run dashboard/app.py
```

Opens at http://localhost:8501.

**Run from the repo root** so Streamlit picks up the `.streamlit/config.toml`
(which enables static-file serving for the row-level brief links).

## Data snapshot

`dashboard/data/repurposing_hypotheses.parquet` is a frozen snapshot of the
full-screen output (`output_full/repurposing_hypotheses.parquet`). To refresh
after re-running the pipeline:

```bash
cp output_full/repurposing_hypotheses.parquet \
   dashboard/data/repurposing_hypotheses.parquet
```

## Deploy to Streamlit Community Cloud

1. Push to a public GitHub repo (the data snapshot is ~8 MB, well under the
   100 MB file limit).
2. https://share.streamlit.io -> "Create app" -> point at the repo.
3. Set the main file to `dashboard/app.py` and Python version to 3.11.
4. Streamlit picks up `dashboard/requirements.txt` automatically.

The app is read-only and has no secrets, so the default public-share settings
are appropriate. For private link-protected hosting, use the Streamlit
"viewer access" controls or deploy behind a reverse proxy.
