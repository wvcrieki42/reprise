# REPRISE dashboard

Public Streamlit companion to the REPRISE drug-repurposing screen. Browses
the 176,272 ranked hypotheses interactively: filter by drug / disease / target,
sort by opportunity, click a row to see the full mechanism / IP / market brief.

## Run locally

```bash
pip install -r dashboard/requirements.txt
streamlit run dashboard/app.py
```

Opens at http://localhost:8501.

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
