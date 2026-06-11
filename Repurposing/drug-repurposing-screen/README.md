# Mechanism-Driven Drug-Repurposing Screen

Screen **every approved drug (US + EU)** against **every disease** and rank all
drug–disease pairs by how strongly the drug's known mechanism of action implicates
a disease it is **not yet indicated for**.

The reasoning chain is the same one demonstrated by hand on auranofin:

```
approved drug ─▶ known targets (MoA) ─▶ STRING functional network
              ─▶ target→disease associations ─▶ subtract known indications
              ─▶ opportunity score ─▶ ranked novel hypotheses
```

This repo is the **engine**. It ships with a small bundled sample dataset so it runs
in seconds with no downloads, and adapters/scripts to run it at full scale against the
standard open databases (ChEMBL, Open Targets, EFO, STRING).

---

## Quick start (demo, no downloads)

```bash
pip install -r requirements.txt
make demo        # or: PYTHONPATH=src python -m repurpose.cli --config config.yaml
make test        # end-to-end smoke test
```

`make demo` writes `output/repurposing_hypotheses.csv` and prints the top hypotheses.
On the bundled data it reproduces recognisable real-world repurposing leads, e.g.
Pioglitazone→NASH, Celecoxib→colorectal cancer, Disulfiram→glioma, Fingolimod→ulcerative
colitis, Imatinib/Sildenafil/Auranofin→IPF — while correctly **excluding** each drug's
existing indications (Auranofin→RA, Sildenafil→PAH, Metformin→T2D … are filtered out).
### Environment setup (macOS)

This project expects a working Python 3 installation and a `python` executable on `PATH`
(the `Makefile` uses `python`, not `python3`).

```bash
# verify runtime
python --version
python3 --version

# install dependencies
python -m pip install -r requirements.txt
```

If `python` is missing but `python3` exists, add a symlink:

```bash
ln -s /usr/local/bin/python3 /usr/local/bin/python
```

If you see NumPy import errors mentioning an incompatible architecture (for example
`x86_64` vs `arm64`), reinstall binary scientific packages for the active runtime:

```bash
python -m pip install --upgrade --force-reinstall --no-cache-dir numpy pandas
```

---

## How it works

| Step | Module | What it does |
|------|--------|--------------|
| 1. Universe | `steps.build_universe` | Keep launched (max_phase 4) drugs approved in US and/or EU. |
| 2. Targets | `steps.build_drug_targets` | Attach each drug's curated MoA targets; optionally expand with STRING functional partners (weighted lower than direct targets). |
| 3. Propagation | `steps.propagate_disease` | For each (drug, disease), combine its targets' disease-association scores into a **mechanistic support** value (noisy-OR / max / sum). |
| 4. Novelty | `steps.add_novelty` | Mark a pair as *known* if the disease is an existing indication or within `ontology_radius` EFO hops of one; otherwise *novel*. |
| 5. Direction | `steps.add_direction` | Combine the drug's `action_type` with the target's disease direction-of-effect into a `direction_factor` (aligned / opposed / unknown). |
| 6. Tissue | `steps.add_tissue` | Gate by whether a driving target is actually expressed in the disease's tissue → `tissue_factor` (expressed / low / absent / unknown). |
| 7. Score | `steps.score` | `opportunity = mechanistic_support^w · novelty^w · direction_factor · tissue_factor`, with an optional drug-level promiscuity penalty; rank descending. |

### Scoring math

For drug *d* and disease *s* with contributing targets *t*:

```
contrib(t, s) = target_weight(t) · assoc(t, s)        # weight 1.0 direct, <1 for network neighbours
mechanistic_support(d, s) = 1 − Π_t (1 − contrib(t, s))     # noisy-OR (default)
novelty(d, s)  = 0            if s is a known indication of d
               = hop/(R+1)    if s is within R ontology hops of one (soft mode)
               = 1            otherwise
opportunity(d, s) = mechanistic_support^w_mech · novelty^w_nov · direction_factor · tissue_factor
                    / sqrt(n_direct_targets(d))      # promiscuity penalty (optional)
```

The promiscuity penalty is **per drug** (uniform across that drug's hypotheses), so it
adjusts comparisons *between* drugs without distorting the ranking of diseases *within* a
drug.

### Directionality model

A target–disease association says the gene *matters*; it does not say whether the drug
moves the target the **therapeutically useful** way. The direction model closes that gap:

```
drug_dir(d, t)       = +1 if action_type raises target activity (AGONIST/ACTIVATOR…)
                       −1 if it lowers it (INHIBITOR/ANTAGONIST/BLOCKER…)        [ChEMBL]
therapeutic_dir(t,s) = +1 if raising activity is therapeutic, −1 if lowering is  [Open Targets genetics]
alignment(d,t,s)     = drug_dir · therapeutic_dir            # +1 aligned, −1 opposed

alignment(d,s)       = Σ_t contrib · alignment(d,t,s) / Σ_t contrib     # over DIRECT targets only
direction_factor     = opposed_factor + (aligned_factor − opposed_factor)·(alignment+1)/2
                     = unknown_factor   when no directional evidence exists
```

`therapeutic_dir` is derived from Open Targets genetic direction-of-effect: e.g. if
loss-of-function variants are *protective*, lowering the target is therapeutic, so an
**inhibitor is aligned** and an **agonist is opposed**. Only direct targets contribute
(a drug's `action_type` says nothing about its network neighbours). The same machinery
flags that auranofin→IPF is *aligned* but auranofin→Friedreich's ataxia could be *opposed*.

### Tissue-expression filter

A mechanism only matters if the drug's target is actually present where the disease is.
For each hypothesis the screen takes the drug's supporting targets, looks up their baseline
expression in the disease's relevant tissue(s), and keeps the best:

```
tissue_expr(t, s)  = max over tissues T relevant to s of  relevance(s,T) · expression(t,T)
tissue_score(d, s) = max over supporting targets t of tissue_expr(t, s)
tissue_factor      = 1.0                         if tissue_score ≥ min_expression  (expressed)
                     scaled toward absent_factor if 0 < tissue_score < min          (low)
                     absent_factor               if measured but absent in tissue   (absent)
                     unknown_factor              if no expression/tissue data       (unknown)
```

On the sample this lifts lung/colon/liver-expressed hypotheses and discounts e.g.
Simvastatin→Alzheimer's, because HMGCR is barely expressed in brain.

Everything is configurable in `config.yaml` (aggregation, thresholds, network on/off,
weights, penalty, and the aligned/opposed/unknown factors). Set `direction.enabled: false`
for a neutral run.

### Engines

Identical results, two backends (`engine: pandas | duckdb`):

- **pandas** — loads everything in memory; ideal up to a few GB.
- **duckdb** — streams the large target–disease table from disk (CSV/Parquet) and does the
  join, noisy-OR aggregation, direction, novelty and scoring in one SQL pass; use it for the
  full Open Targets association set. The smoke test asserts the two engines agree exactly.

---

## Running at full scale (all approved drugs)

```bash
bash scripts/download_data.sh        # ChEMBL + Open Targets + EFO (+ optional STRING/DrugBank)
python scripts/build_full_tables.py  # materialise data/full/*.csv in canonical schema
# edit config.yaml: set `mode: full` and repoint paths to data/full/*
make full
```

**Data sources (all open):**

| Table | Source | Notes |
|-------|--------|-------|
| Approved drugs + region | ChEMBL `molecule_dictionary` (max_phase 4); merge FDA Orange Book / EMA list for accurate EU flag | ~4–5k launched drugs |
| Drug→target (MoA) | ChEMBL `drug_mechanism` | curated mechanism targets |
| Drug→indication | ChEMBL `drug_indication` (EFO-mapped) | the novelty subtraction set |
| Target→disease | Open Targets `associationByOverallDirect` | genetics + expression + literature evidence, score 0–1 |
| Target direction | Open Targets `evidence` (variant effect × direction-on-trait) | therapeutic direction-of-effect, `±1` |
| Target expression | Open Targets `baselineExpression` (GTEx/HPA) | per-tissue expression for the tissue gate |
| Disease→tissue | curated / EFO–UBERON cross-refs | maps each disease to its relevant tissue(s) |
| Gene names | Open Targets `targets` (`approvedName`) | lead-target full name for the output |
| Disease ontology | EFO `.obo` | parent/child hops for the novelty radius |
| Network | STRING v12 (API or full human edge file) | functional neighbourhood expansion |

The adapters in `src/repurpose/sources/adapters.py` contain the exact SQL/parsers for each
(including `opentargets_target_direction` for the directionality table).

**Scale note.** The full join is `drugs × targets × diseases` ≈ tens of millions of
candidate pairs but extremely sparse after thresholding. The pandas engine handles it in a
few GB of RAM; for the full Open Targets association set, set `engine: duckdb` to stream the
big table from disk and run the join/aggregation/scoring in SQL — identical results,
out-of-core.

---

## Output

`output/repurposing_hypotheses.csv` (and `.parquet`), one row per ranked hypothesis:

| column | meaning |
|--------|---------|
| `rank`, `opportunity` | final priority |
| `drug_id`, `drug_name`, `modality` | the drug |
| `efo_id`, `disease_name` | the candidate indication |
| `lead_target` | top-contributing target (gene symbol) |
| `lead_target_name` | full target/protein name |
| `lead_target_genecards` | GeneCards URL for the lead target |
| `mechanistic_support` | strength of the MoA→disease link (0–1) |
| `novelty`, `novelty_status` | novel / known_exact / known_related_hN |
| `direction_factor`, `direction_status` | aligned / opposed / mixed / unknown (therapeutic-direction gate) |
| `tissue_factor`, `tissue_status`, `tissue_evidence` | expressed / low / absent / unknown, and the tissue |
| `n_targets`, `n_drug_targets` | supporting targets for this disease / total direct targets of the drug |
| `evidence_targets` | the targets (and contributions) that drove the hypothesis |
| `data_version` | source versions (ChEMBL / Open Targets / STRING / EFO) for provenance |

A sidecar `output/run_metadata.json` records the run date, engine, mode, data versions,
thresholds and hypothesis count.

---

## Repo layout

```
config.yaml                  run configuration (demo by default)
Makefile                     install / demo / test / full
requirements.txt
src/repurpose/
  config.py                  typed config loader
  ontology.py                EFO parent/child navigation (novelty radius)
  steps.py                   the pipeline steps (pure functions)
  pipeline.py                orchestration + engine dispatch
  cli.py                     `python -m repurpose.cli`
  sources/
    loaders.py               canonical-table readers (CSV/Parquet)
    adapters.py              ChEMBL / Open Targets / EFO extractors (full mode)
    string_api.py            STRING network expansion (file or live API)
  backends/
    duckdb_engine.py         out-of-core SQL engine (engine: duckdb)
data/sample/                 bundled demo dataset (illustrative, not authoritative)
scripts/
  download_data.sh           fetch the full bulk datasets
  build_full_tables.py       convert dumps -> canonical tables
tests/test_smoke.py          end-to-end test + logic assertions
```

---

## Limitations (read before believing any single hypothesis)

- **Association ≠ direction — partially handled.** A target–disease association tells you the
  gene *matters*, not whether **inhibiting** vs activating it helps. The directionality model
  addresses this where genetic direction-of-effect exists, but coverage is incomplete: many
  pairs fall to `unknown` and are only mildly discounted, not vetoed. Where evidence is thin or
  the disease mechanism is non-monotonic (e.g. auranofin's redox effects in Friedreich's ataxia),
  the gate can still be wrong. **Every top hypothesis needs a human mechanistic sanity check.**
- **Target sets are incomplete** and ignore off-targets, active metabolites, and tissue
  exposure / pharmacokinetics. A mechanistic match still has to reach the diseased tissue at a
  tolerated dose.
- **The bundled `data/sample` is illustrative**, hand-built to exercise the code — do not cite
  its numbers. Real conclusions require the full ChEMBL/Open Targets run.
- **Output is hypothesis-generating.** It ranks where to look; it does not establish efficacy.
  Confirm with cell and animal models before any clinical claim.

## License / data terms
Code: use freely. Respect the individual licences of ChEMBL (CC BY-SA), Open Targets
(CC0), EFO, STRING, and DrugBank (separate licence) when using their data.
```
