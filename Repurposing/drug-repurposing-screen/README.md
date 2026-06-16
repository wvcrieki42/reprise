# Drug Repurposing Screen
Mechanism-driven screen to rank potential novel drug–disease pairs from known drug targets, target–disease evidence, novelty filtering, and configurable scoring.

**Public dashboard: https://reprise.streamlit.app** — browse the 176,272-row ranked output interactively, no install required.

## Project structure
```text
drug-repurposing-screen/
├── Makefile
├── README.md
├── requirements.txt
├── config.yaml
├── config.full.yaml
├── src/
│   └── repurpose/
│       ├── cli.py
│       ├── config.py
│       ├── ontology.py
│       ├── pipeline.py
│       ├── steps.py
│       ├── backends/
│       │   └── duckdb_engine.py
│       └── sources/
│           ├── adapters.py
│           ├── loaders.py
│           └── string_api.py
├── data/
│   └── sample/
├── scripts/
│   ├── build_full_tables.py
│   └── download_data.sh
└── tests/
    └── test_smoke.py
```

## Prerequisites
- Python 3.11+ installed.
- `python` and `pip` available on `PATH`.

Check:
```bash
python --version
python3 --version
python -m pip --version
```

If `python` is missing but `python3` exists:
```bash
ln -s /usr/local/bin/python3 /usr/local/bin/python
```

## Setup
Install dependencies:
```bash
python -m pip install -r requirements.txt
```

## Run the project
Demo run:
```bash
make demo
```

Smoke tests:
```bash
make test
```

Full run (after downloading/building full datasets):
```bash
bash scripts/download_data.sh
python scripts/build_full_tables.py
make full
```

## Output
- `output/repurposing_hypotheses.csv`
- `output/run_metadata.json`

## Web dashboard
Browse the ranked hypotheses interactively:
```bash
pip install -r dashboard/requirements.txt
streamlit run dashboard/app.py
```
See `dashboard/README.md` for deployment to Streamlit Community Cloud.

## Latest demo verification
Most recent local validation run:

```bash
make demo
```

Result:
- Exit code: `0` (success)
- Universe size: `12` approved US/EU drugs
- Final ranked hypotheses: `35`
- Artifacts written:
  - `output/repurposing_hypotheses.csv`
  - `output/run_metadata.json`

## Common troubleshooting
### NumPy architecture mismatch (Apple Silicon)
If you get import errors mentioning incompatible architecture (`x86_64` vs `arm64`):
```bash
python -m pip install --upgrade --force-reinstall --no-cache-dir numpy pandas
```

### Missing `python` command
The Makefile uses `python`. Ensure `python` resolves in shell:
```bash
command -v python
```