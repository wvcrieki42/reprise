# A Mechanism-Driven Drug-Repurposing Screen with Integrated Evidence Layers, Intellectual-Property Signals, and Operational Outputs

*Wim Van Criekinge, BioBix Group, Department of Mathematical Modelling, Statistics and Bioinformatics, Ghent University, Belgium*

---

## Abstract

**Background.** Computational drug-repurposing screens often surface candidates that look plausible in isolation but fail at the next layer of scrutiny: an indication may already be on-label, the underlying biology may be dose-impractical, the intellectual property may have expired, or no clinician-leader is available to drive a trial. We present an end-to-end screen that scores approved drugs against the full disease ontology and layers on enrichment passes — directionality, tissue expression, orthologous-gene model-organism evidence, Reactome pathway co-membership, FDA Orange Book intellectual-property runway, Orphanet rare-disease prevalence, PubMed/Europe PMC/ClinicalTrials.gov investigation prior, and an automated identification of regional Key Opinion Leaders — to produce a single ranked output that is interpretable end-to-end and immediately operational.

**Methods.** Mechanistic support per (drug, disease) pair was computed as the noisy-OR over drug-target × target-disease contributions, using ChEMBL approved drugs and Open Targets (24.06) target-disease association scores. Novelty was inferred from the EFO ontology relative to known ChEMBL drug-indication assignments, with a target-class indication rollup applied to close formulation-coverage gaps. Direction-of-effect was inferred from Open Targets genetic evidence after filtering to human-genetics datasources only (IMPC mouse-knockout signals were excluded after we showed they bias monogenic-disorder direction calls). Phylogenetic and pathway-level evidence were added as multiplicative boosts. A two-pass literature enrichment hits PubMed/Europe PMC/ClinicalTrials.gov and an optional Lens.org patent backend; counts were combined into a log-saturating investigation_prior that damps well-trodden hypotheses. Formulation variants were collapsed up the ChEMBL active-ingredient hierarchy. Severity flag damps the receptor-LoF-plus-agonist class. Each top hypothesis was matched to one US and one EU Key Opinion Leader from PubMed authorship, with affiliation-based regional preference (Boston / Cambridge-MA / Bay Area for US; Belgium / Netherlands / France-Germany / UK for EU) and h-index resolved via Semantic Scholar. The pipeline emits a 36-column CSV plus one-page PDF deal memos.

**Results.** On Open Targets release 24.06 the full screen produced 176,265 ranked hypotheses in approximately 90 seconds of compute (DuckDB engine, M-series Mac), with the literature pass adding ~10 minutes for the top 5,000 hypotheses when an NCBI API key was supplied. Backtest validation against 19 curated known repurposing successes (sildenafil/PAH, thalidomide/multiple myeloma, raloxifene/osteoporosis, methotrexate/RA, etc.) yielded an 84% hit rate at mech_support ≥ 0.3, with all three misses traceable to specific Open Targets data-coverage gaps for rare or peripheral indications. Severity damping pushed all 32 receptor-LoF + agonist trap hypotheses out of the top 100. Per-substance grouping reduced the output from 200,000 raw drug-disease rows to 175,759 substance-disease rows, with formulation variants reported in a `variant_names` column. The Orphanet bridge added structured prevalence data for ~5,100 rare diseases, lifting market-data coverage from 18,229 to 32,292 hypotheses (a 77% increase) and orphan-flagging from 5,607 to 19,355 (a 3.5-fold lift). The Orange Book join populated patent / exclusivity / generic-status columns for ~62,000 hypotheses — biologics correctly remained NA, in line with their Purple Book classification.

**Conclusion.** A repurposing screen is most useful when its mechanistic prediction is paired with the regulatory, economic, and human-network signals needed to actually advance a candidate. The pipeline presented here demonstrates that these layers can be integrated end-to-end on free public data, with measured precision against known successes, and with operational outputs (PDF deal memos with identified KOLs) suitable for immediate clinical-partnership outreach.

---

## 1 Introduction

Pharmacological repurposing — bringing an already-approved drug to a new indication — has produced some of the most consequential changes in modern medicine: sildenafil from angina to pulmonary arterial hypertension, thalidomide from morning sickness to multiple myeloma, methotrexate from oncology to rheumatoid arthritis, propranolol from hypertension to infantile hemangioma [refs]. The economic and clinical case for systematic repurposing is now well-established: an approved drug already carries through safety and ADME validation, shortens development time by years, and reduces capital-at-risk by an order of magnitude relative to de-novo drug discovery [refs].

In silico repurposing screens have proliferated alongside the public biomedical databases that make them feasible. PREDICT used drug-drug and disease-disease similarity matrices to suggest novel pairings [Gottlieb 2011]. DrugCentral aggregates regulatory and mechanistic curation [Avram 2021]. The Connectivity Map and L1000 platforms exploit transcriptional signatures to relate drugs and diseases [Subramanian 2017]. More recently, network-pharmacology approaches and graph neural networks have extended the scope toward heterogeneous evidence integration [refs].

A persistent weakness across this family of tools, however, is that they typically end at a mechanism score. Useful as that is, a mechanism score alone cannot answer the questions that determine whether a candidate advances:

- Is the indication actually novel relative to ChEMBL annotation, or already on-label?
- Does the drug push the target in the direction therapeutically appropriate for the disease?
- Is the receptor actually expressed in the disease tissue?
- For monogenic disorders, is the drug an agonist for a target the disease has already broken?
- What is the substance's intellectual-property runway, and is there any generic competition?
- What is the addressable patient population, and does the indication qualify for orphan exclusivity?
- Has the literature already saturated the hypothesis, or is it open ground?
- Which clinical-researcher is best positioned to lead a trial?

The work presented here is an attempt to answer all of these in one integrated pipeline, with measured precision against known repurposing wins and with operational outputs (per-hypothesis PDF deal memos including identified Key Opinion Leaders) ready for clinical-partnership outreach. We use only free public data: Open Targets [Carvalho-Silva 2019], ChEMBL [Mendez 2019], Reactome [Gillespie 2022], STRING [Szklarczyk 2023], EFO [Malone 2010], FDA Orange Book, Orphanet [Pavan 2017], PubMed E-utilities, Europe PMC, ClinicalTrials.gov v2 API, and Semantic Scholar [Wang 2020].

The pipeline is open source and reproducible from a documented `download_data.sh` step; one full run produces ~176,000 ranked hypotheses in ~90 seconds of compute. We measure its precision by backtesting against 19 curated known repurposing successes, achieving an 84% hit rate (16/19) at a mechanism-support threshold of 0.3.

The remainder of this manuscript describes the methods, results, and limitations in detail.

---

## 2 Methods

### 2.1 Data sources

The pipeline draws from the following free public sources, all integrated at fixed release versions:

- **ChEMBL 37** [Mendez 2019]: approved drugs, their direct molecular targets via mechanism-of-action annotations, drug indication mappings, and the molecule_hierarchy table for active-ingredient rollup. Accessed via the SQLite distribution.
- **Open Targets 24.06** [Carvalho-Silva 2019]: target-disease association scores (`associationByOverallDirect`), disease ontology metadata (`diseases`, including synonyms and dbXRefs), gene-name and ensembl-id maps from the `targets` parquet, baseline expression by tissue (`baselineExpression`), and the full evidence parquet (`evidence`) for direction-of-effect inference, animal-model phylogenetic evidence, and Reactome pathway memberships per gene.
- **Reactome** (via Open Targets' integrated `targets.pathways` field): gene-pathway memberships used for indirect-mechanism evidence.
- **STRING v12.0** [Szklarczyk 2023]: protein-protein interaction network, downloaded as protein.links and protein.info, then converted to a symbol-keyed edge list (~1.9M edges at min_score ≥ 400) for optional neighbour expansion of drug-target sets.
- **EFO 2024-06** [Malone 2010]: disease ontology graph used to expand "known indication" via parent-child relationships up to a configurable hop radius.
- **FDA Orange Book** (May 2026 monthly release): patent and exclusivity expirations and ANDA (generic) application counts per active ingredient, parsed from the three tilde-separated text files (products / patent / exclusivity).
- **Orphanet** (en_product9_prev.xml, CC-BY-4.0): structured prevalence categories for ~6,400 rare diseases, bridged to Open Targets MONDO/EFO identifiers via the OT `dbXRefs` field.
- **PubMed E-utilities, Europe PMC REST, ClinicalTrials.gov v2 API**: literature counts per (target, disease) query, with synonym expansion, used for the investigation_prior signal.
- **Lens.org Patent API**: optional patent-search backend.
- **Semantic Scholar Author API** [Wang 2020]: h-index resolution for identified Key Opinion Leaders.

A single curated CSV (`disease_prevalence.csv`) supplies US prevalence figures for 122 major diseases hand-sourced from CDC, SEER, foundation estimates, and similar authoritative sources. The Orphanet automated extraction supplements this with ~5,100 additional rare-disease rows; on EFO collision the curated entry wins.

### 2.2 Mechanistic support

For each (drug, disease) pair, mechanistic support was computed as in earlier mechanism-driven screens [refs]:

1. Drug-target edges were taken from ChEMBL mechanism-of-action annotations. Each direct target carried `target_weight = 1.0`. If `network.enabled`, the target set was expanded via the STRING PPI graph (default `min_confidence = 0.7`, `max_partners = 10`); expanded neighbours carried `target_weight = neighbor_weight × STRING_score` with `neighbor_weight = 0.5` by default.
2. For each (drug, target, disease) edge, the per-edge contribution was

```
contrib = clip(target_weight, 0, 1) × clip(assoc_score, 0, 1)
```

where `assoc_score` is the OT association score, filtered by `min_assoc ≥ 0.1`.

3. Per (drug, disease), edge contributions were aggregated by a configurable function — by default, the **noisy-OR**:

```
mech_support = 1 − ∏ᵢ (1 − contribᵢ)
```

with `max` and capped-sum aggregators available as alternatives.

This is the same per-pair score used by direction-aware mechanism screens dating back to PREDICT-style approaches but with the noisy-OR specifically chosen because it captures "any one of these is sufficient" semantics that match the intuition of mechanistic-redundancy in the disease's biology.

### 2.3 Novelty

Novelty per (drug, disease) was defined as ontological distance from the drug's nearest known indication. The pipeline takes ChEMBL `drug_indications` as the baseline known set (after a colon→underscore normalisation that we found necessary, because ChEMBL exports MONDO and HP-coded indications with colons while Open Targets uses underscores). Around each known EFO, the EFO ontology graph was walked outward by a configurable radius (default 1), with the novelty score assigned by

```
hop 0:  novelty = 0.0, status = known_exact
hop h > 0, soft = true:  novelty = h / (radius + 1), status = known_related_h{h}
not reached:  novelty = 1.0, status = novel
```

A subtle but important refinement is the **target-class indication rollup**. ChEMBL annotates indications on individual drug rows, but biologics often have one row per formulation (e.g., `INSULIN ASPART` vs. `INSULIN ASPART PROTAMINE RECOMBINANT`) and indications are patchy across them; the same is true for selective vs. broad members of small-molecule classes (e.g., `SILODOSIN` vs. `DOXAZOSIN`, both ADRA1 antagonists). Without rollup, the silodosin/hypertension pairing surfaced spuriously at the top of the screen. The rollup rule we adopted is: drug A inherits drug B's indications iff A's frozenset of (direct target_symbol, action_type) pairs is a SUBSET of B's. Selective drugs inherit from broader drugs; the reverse is not asserted. STRING neighbours are excluded from class keys to prevent inflation.

### 2.4 Directionality

Direction-of-effect per (target, disease) was inferred from Open Targets genetic evidence following the formulation in the OT documentation: variant effect (LoF vs. GoF) combined with direction-on-trait (risk vs. protective) implies a therapeutic direction in {+1, −1} on whether the target activity should be increased or decreased. We found during development that the IMPC mouse-knockout source heavily biases this calculation for monogenic disorders — for example, for EPOR + primary familial polycythemia (a gain-of-function truncating disorder), IMPC contributed 16 LoF-risk rows producing a misleading "activate EPOR" inference. After excluding IMPC and other non-human-genetics sources (cancer_biomarkers, chembl, europepmc, etc.) and keeping only `ot_genetics_portal`, `gene_burden`, `eva`, `gene2phenotype`, `clingen`, `genomics_england`, and `orphanet`, the direction table contracted from 1,019,063 rows to 113,339 rows and the EPO-polycythemia false positive dropped 275 ranks.

The per (drug, disease) direction factor was computed as the contribution-weighted alignment of (drug action direction) × (therapeutic direction) over informative direct-target edges, mapped to `[opposed_factor, aligned_factor]` via a configurable linear interpolation. The default opposed factor was 0.15, aligned 1.0, and unknown (no informative edges) 0.6. The score is multiplicative.

### 2.5 Tissue gate

Tissue evidence asked whether a driving target is actually expressed where the disease manifests. From OT baseline expression and a curated disease → tissue map (built by mining OT therapeutic-area metadata plus a small set of disease-name overrides), per-target tissue expression scores were maximised across the disease's relevant tissues. Per (drug, disease), the best target's tissue score determined whether the hypothesis was classified `expressed`, `low`, `absent` (measured zero), or `unknown` (no measurement), with corresponding factors (1.0, ramp, `absent_factor=0.3`, `unknown_factor=0.7`). The factor is multiplicative.

### 2.6 Phylogenetic boost

For each disease's strong-association targets, OT model-organism evidence (`datatypeId == animal_model`, defaulting to IMPC mouse knockouts via PhenoDigm) was aggregated as a max-over-rows score in `[0, 1]`. The per (drug, disease) phylo_factor was computed as

```
phylo_factor = 1 + boost_factor × max(phylo_score over drug's targets)
```

with default `boost_factor = 0.5` (so the maximum boost is 1.5×). This is asymmetric by design: present evidence boosts, absent evidence does not penalise (the corresponding mouse knockout may simply not have been run).

### 2.7 Pathway evidence

Indirect-mechanism evidence was added via Reactome pathway co-membership, as exposed in the OT `targets.pathways` field. For each (drug, disease) pair, we counted the unique Reactome pathways (filtered to those with ≤ 80 genes, to exclude generic supercategories) that contained both a direct target of the drug AND a strong-association target of the disease. Crucially, we required the drug-side target and disease-side target to **differ** (`indirect_only = true`), so the pathway pass only credits cases the direct-target signal misses; this prevents double-crediting of (drug, disease) pairs where the drug already directly hits a disease-driving target. The per (drug, disease) pathway_factor was

```
pathway_score = log1p(n_overlap) / log1p(saturation),    saturation = 5
pathway_factor = 1 + boost_factor × pathway_score,        boost_factor = 0.3
```

again asymmetric (absence does not penalise).

### 2.8 Final opportunity score and substance grouping

The opportunity score, by default, is

```
opportunity = mech^w_mech × novelty^w_novelty × direction_factor × tissue_factor × phylo_factor × pathway_factor / sqrt(n_drug_targets)
```

with weights `w_mech = w_novelty = 1.0` and a promiscuity penalty of `sqrt(n_drug_targets)` (default on). After scoring, the substance-grouping step collapses formulation variants up to their ChEMBL active-ingredient parent via the `molecule_hierarchy` table — `PIOGLITAZONE HYDROCHLORIDE → PIOGLITAZONE`, `ROSIGLITAZONE MALEATE → ROSIGLITAZONE`, etc. — keeping the highest-opportunity row per (substance, disease) as representative and listing the variant names in a `variant_names` column.

### 2.9 Severity flag

A narrow heuristic catches the "give the agonist to the broken receptor" trap: when (a) the `disease_gene_match` flag indicates that a gene named in the disease is also a direct target of the drug, (b) the disease name contains a severity keyword (`deficiency`, `complete absence`, `Donohue`, `Rabson-Mendenhall`, etc.), and (c) the drug's action_type on the named target is in the agonist class, the hypothesis is flagged `severe_loF_agonist` and its opportunity is multiplied by `(1 − damp_factor)` with `damp_factor = 0.7` by default. The heuristic is deliberately conservative: TZDs for PPARG-related familial partial lipodystrophy and calcimimetics for familial hypocalciuric hypercalcemia 1 are NOT flagged (their disease names lack severity keywords), preserving these as legitimate hypotheses where partial residual function is rescuable by agonist therapy.

### 2.10 Literature investigation prior

For the top-N hypotheses (default 5,000), we query PubMed, Europe PMC, ClinicalTrials.gov v2, and optionally Lens.org's patent API. Queries are constructed as

```
(target term OR target gene name OR synonyms) AND (disease name OR exact synonyms)
```

with synonyms pulled from OT `targets.approvedName` and `diseases.synonyms.hasExactSynonym` (broad/related/narrow synonyms intentionally excluded — they lower precision). Per source, the hit count was normalised via log-saturation:

```
norm = min(log1p(count) / log1p(saturation), 1.0)
```

and the per-source scores combined via a weighted mean to give `investigation_prior ∈ [0, 1]`. The trial signal was weighted 2× the literature signals. When `scoring.w_investigation > 0`, opportunity was damped by `(1 − w × investigation_prior)` and the output re-ranked. By default `w_investigation = 0.5`, which damped well-trodden combinations enough that the top of the screen surfaced genuinely under-explored hypotheses while preserving the original mechanism-driven scoring for the inspection-only case.

### 2.11 Market sizing

Per disease, US patient counts and prevalence per 100,000 were sourced from two layered datasets: a manually curated CSV of 122 major US diseases (with provenance per row from CDC, SEER, NIDDK, and patient-foundation estimates), supplemented by a script that parsed Orphanet's bulk prevalence XML for the rare-disease long tail (~5,100 disorders). On EFO-id collision the curated entry wins. Orpha codes were bridged to OT's primary disease IDs (mostly `MONDO_*`) via the OT `diseases.dbXRefs` field, so the rare-disease prevalence rows key by the same id family the rest of the pipeline uses. The `is_orphan` boolean flag was derived against the FDA Orphan Drug Act 200,000-US-patient threshold and propagated as a nullable Boolean (NA-when-unknown, not False, by design).

### 2.12 Intellectual-property signals

The FDA Orange Book was parsed to produce per-active-ingredient rows containing latest patent year, latest exclusivity year, the maximum (LOE year), and a Boolean for generic availability (any ANDA filing). Combination products were split on the semicolon delimiter and attributed to each component. Salt suffixes were stripped (`PIOGLITAZONE HYDROCHLORIDE → PIOGLITAZONE`) for matching against the substance_name post-collapse. Biologics are outside the Orange Book scope (they live in the Purple Book) and consequently return NA — a correct outcome given the question "has this gone generic / how long until LOE."

### 2.13 Key Opinion Leader identification

For the top-N hypotheses with a literature signal (default 100), we ran a separate PubMed esearch + efetch over the (target × disease) query, parsed `AuthorList` and `AffiliationInfo` fields, classified affiliations by region using a curated keyword set (with explicit `Cambridge, MA` vs. `Cambridge UK` disambiguation), scored affiliations by city / institution preference, and picked the highest-scoring author per region. Affiliation preferences: for the US, Boston / Cambridge MA / Bay Area (San Francisco, Stanford, Berkeley, Palo Alto, UCSF) > other US; for Europe, Belgium > Netherlands > France / Germany / Luxembourg > United Kingdom > other Europe. Reflecting the host institution of this work, the EU preference list was tuned to identify candidates close to Ghent University.

h-index resolution used Semantic Scholar's free Author API. A common author-name convention in Semantic Scholar uses abbreviated first names (`E. Topol` vs. PubMed's `Eric Topol`), which a naïve `==` match misses. We required the surname to match exactly and the first-name first-letter to agree; among candidates passing this filter, we selected the highest h-index — when multiple authors share a surname plus an initial, the most-cited author is overwhelmingly the more likely KOL. This sort rule was validated against a case in development where `Boris Keren` was being matched to a junior author with h=1 over the established geneticist with h=56.

### 2.14 Combination-therapy companion finder

For each top-N primary hypothesis, the pipeline searches for the 1–2 most synergistic companion drugs. For each strong disease target the primary doesn't hit, candidate substances hitting that target are evaluated. The combo mech_support is the noisy-OR over the union of primary and candidate contribution edges, and the synergy is `mech_combo - mech_primary`. Candidates with target overlap greater than 50% are rejected as same-class redundancy (we want different mechanisms, not duplicates). The top 1–2 candidates by synergy are reported.

### 2.15 Validation methodology

A curated YAML at `data/curated/repurposing_validation.yaml` lists 19 known successful repurposing cases (sildenafil/PAH, thalidomide/MM, raloxifene/osteoporosis, methotrexate/RA, hydroxychloroquine/SLE, spironolactone/heart-failure, memantine/AD, tofacitinib/UC, naltrexone/alcohol, duloxetine/fibromyalgia, prazosin/PTSD, aripiprazole/MDD, colchicine/pericarditis, plus several others). For each case, the validation script computes mech_support exactly as the pipeline does — noisy-OR over drug-target × target-disease at `min_assoc = 0.1` — using all ChEMBL IDs that share the drug name (parent + salt variants — necessary because ChEMBL's mechanism-of-action coverage is patchy at the parent level). A case is scored HIT at `mech_support ≥ 0.3`. Failures are categorised: DRUG_NOT_IN_UNIVERSE (drug missing from approved list), DISEASE_NOT_IN_OT (disease ID not found), or MISS (mech_support below threshold despite drug and disease being present).

### 2.16 Operational outputs

A PDF deal-memo generator (`scripts/build_deal_memos.py`) renders one-page A4 memos per (top hypothesis, region) combination using ReportLab. Each memo weaves the mechanism signals into prose (mech_support level, novelty status, direction, tissue, phylo, pathway), translates the Orange Book columns into an IP-runway paragraph ("LOE 2037, no generic" or "LOE 2027, has generic"), expresses the market opportunity ("approximately 16,650 US patients, below the 200,000 Orphan Drug Act threshold"), introduces the identified KOL with credentials ("h-index 56, 1 PubMed publication co-mentioning the target and indication"), and concludes with a proposed three-way collaboration model (KOL → bioinformatics analytics partner → IP-holder substance owner under an investigator-initiated trial framework). A profile YAML supplies the affiliation, principal investigator name and credentials, and bioinformatics-pitch bullets used in each memo, with command-line overrides for one-off customisation.

---

## 3 Results

### 3.1 Full-run statistics

A single end-to-end run on Open Targets release 24.06 with all enrichment layers active and the literature pass enabled at `top_n = 5,000` produced:

```
universe (approved US/EU):    3,996 drugs
drug-target edges:           38,649 (10,196 direct + 28,453 PPI neighbour)
target-disease associations: 376,810 (at min_assoc ≥ 0.1)
ranked hypotheses:          200,000 (capped by scoring.top_n)
after substance collapse:   175,759 rows
pathway-boosted:            128,313
phylo-boosted:               13,740
with US market data:         32,292 (curated + Orphanet)
orphan-flagged (<200k):      19,355
with Orange Book IP data:    62,210
severity-damped:                 32 (all pushed below rank 100)
```

DuckDB compute time was approximately 21 seconds for the core join + scoring, plus 6 seconds for pathway, 4 seconds for substance collapse, ~10 minutes for the literature pass (5,000 (target, disease) queries × 3 sources × 10 req/s with an NCBI API key), 1 second for market enrichment, 1 second for Orange Book enrichment, and approximately 1 minute for the KOL pass at `top_n = 100`. Total wall clock with everything active was approximately 12 minutes; with the literature cache fully warm, ~30 seconds.

### 3.2 Backtest precision

Of 19 curated known repurposing successes, 16 (84%) scored HIT at `mech_support ≥ 0.3`. The three failures were:

- **Minoxidil → alopecia**. Minoxidil's targets ABCC9 and KCNJ11 (K-ATP channel opener) are present, and minoxidil's hair-growth effect is a well-known peripheral effect, but the Open Targets target-disease layer does not strongly link K-ATP genes to alopecia as a disease association — consistent with alopecia being characterised at a follicle-biology level OT does not curate.
- **Propranolol → capillary infantile hemangioma**. Propranolol's targets ADRB1 and ADRB2 are in the screen, but OT's target-disease support for adrenergic-receptor associations with infantile hemangioma is below the `min_assoc = 0.1` filter — consistent with hemangioma being a rare pediatric vascular condition with limited GWAS or burden-test signal.
- **Acetazolamide → idiopathic intracranial hypertension**. The disease is absent from OT entirely (DISEASE_NOT_IN_OT), so the score cannot be computed.

All three failures are explainable as Open Targets data-coverage gaps rather than pipeline-logic failures. The screen scores correctly on the cases where OT has the underlying target-disease evidence.

A particularly informative case is **Metformin → polycystic ovary syndrome**, which scored HIT with `mech_support = 1.00`. The literature most commonly attributes metformin's PCOS effect to AMPK activation via PRKAA1/PRKAA2, but the screen surfaced the connection via mitochondrial complex I genes — `MT-ND2, MT-ND4, MT-ND4L, MT-ND5, NDUFA13, NDUFAB1, NDUFB8, NDUFB9` — rather than the AMPK subunits. This is actually defensible biology: metformin's primary action is mitochondrial complex I inhibition, with AMPK activation as a downstream consequence [refs]. The screen surfaced the more proximal mechanism than the literature heuristic and is, by that standard, correct on a finer level than the validation YAML's `expected_mechanism_targets` field anticipated. The validation script flags this kind of case explicitly as a "surprising hit" for follow-up review.

### 3.3 Severity flag impact

Without severity damping, the top 30 of the full run was contaminated by an "insulin → hyperinsulinism due to INSR deficiency" cluster (20 rows, ranks 11–30, all sharing the same INSR direct hit and mech_support 0.73). Conceptually these hypotheses are non-actionable: when INSR is functionally null, exogenous insulin agonist will not engage the receptor. With the severity flag enabled (`damp_factor = 0.7`), all 20 rows received a 70% opportunity penalty and dropped to ranks 37,906–37,926 — effectively out of the screen.

Critically, the heuristic preserved adjacent partial-LoF cases that *are* legitimate repurposing strategies: thiazolidinediones for PPARG-related familial partial lipodystrophy (a real clinical strategy where TZD agonist activity rescues residual receptor function) remained at rank 30, and calcimimetics for familial hypocalciuric hypercalcemia 1 (where partial CaSR function is preserved) remained at ranks 19–20. The discrimination rests on the disease name itself: "deficiency" in INSR-deficiency hyperinsulinism implies a severe-LoF state where agonism is futile; "partial lipodystrophy" or "hypocalciuric" does not.

### 3.4 Top-30 read-out

After the full enrichment stack, the top 30 of the run was a mix of off-label clinical reality and orphan-disease opportunities with substantial IP runway. Selected highlights:

- **Voxelotor → 4 distinct hemoglobin disorders** (Heinz body anemia, hereditary persistence of fetal hemoglobin, hemoglobin M disease, hemoglobin H disease). Voxelotor is approved for sickle cell disease (HbS allosteric stabiliser); the cross-indications surface naturally from the shared HBB/HBA target set. Orange Book LOE 2037, no generic competition, ~3 NDA filings — a strong intellectual-property profile.
- **BRAF inhibitors (encorafenib, dabrafenib, vemurafenib) and KRAS inhibitors (sotorasib, adagrasib) → cardiofaciocutaneous syndrome**. CFC syndrome is a developmental RASopathy driven by BRAF/MEK/KRAS pathway dysregulation. The screen connects the targeted-oncology agents to the developmental RASopathy via the shared signalling pathway. LOE 2032–2043 across the cluster, no generic competition. Pathway evidence contributed substantially.
- **Bromazepam → developmental and epileptic encephalopathy**. Bromazepam is a benzodiazepine; DEE is a class of severe pediatric epilepsies frequently driven by mutations in GABA-A receptor subunits (GABRA1, GABRB3, GABRG2). The mechanism is direct (benzo → GABA-A enhancement). Off-label use of benzodiazepines for refractory pediatric epilepsy is reported clinically.
- **Follitropin → 46,XX gonadal dysgenesis**. FSH for ovarian dysgenesis — direct mechanistic plausibility through FSHR signalling, with practical caveats about partial vs. complete dysgenesis.
- **Afamelanotide → oculocutaneous albinism type 6**. Afamelanotide is an α-MSH analogue; oculocutaneous albinism is a tyrosine-pathway pigmentation disorder. LOE 2033, no generic. Rank 9, ~166 US patients (ultra-rare orphan).

These hypotheses are not new to specialists in their respective fields; what is novel is that they emerge from a single mechanism-driven pipeline over the full disease ontology, with the same scoring layer that surfaces the known repurposing successes in the backtest set.

### 3.5 Combination-therapy companions

At `top_n = 200`, the combination-therapy companion finder identified at least one mechanistically complementary partner for 61 primary hypotheses (synergy ≥ 0.1, target overlap ≤ 0.5). Several of the surfaced combinations reproduce established or actively-investigated clinical strategies:

- **Dabrafenib + Binimetinib (and Encorafenib + Binimetinib, Vemurafenib + Binimetinib) → cardiofaciocutaneous syndrome**. The BRAF + MEK doublet recapitulates the standard-of-care combination in BRAF V600E metastatic melanoma, surfaced here from first principles in the developmental RASopathy context. Synergy 0.15. Bridge target: MAP2K2.
- **Sotorasib / Adagrasib + Binimetinib → cardiofaciocutaneous syndrome**. The KRAS + MEK combination, mirroring the BRAF analogue, with synergy 0.17.
- **Angiotensin II + Aliskiren → renal tubular dysgenesis**. Angiotensin agonism + renin inhibition — both push on the same renin-angiotensin axis that drives the developmental defect, in opposite directions. Bridge target: REN. Synergy 0.21.
- **Bromazepam + Orphenadrine → developmental and epileptic encephalopathy**. Benzodiazepine GABA-A potentiation + NMDA antagonism, two different inhibitory-signalling strategies for refractory epilepsy. Synergy 0.48 — the highest in the top 30.
- **Pinacidil / Minoxidil + Vernakalant → familial atrial fibrillation**. K-ATP channel opener (cardiac repolarisation) + atrial-specific K⁺ blocker. Synergy 0.26. Bridge target: KCNJ5.
- **Follitropin + Esterified estrogens → 46,XX gonadal dysgenesis**. FSHR activation + ESR2 — both contributing to follicular development. Synergy 0.11.
- **Setmelanotide / Bremelanotide + Metformin → type 2 diabetes mellitus**. MC4R agonist (anti-obesity / energy expenditure) + metformin (complex I inhibition / glycemic). Synergy 0.49. The NDUFAB1 bridge identifies metformin's mitochondrial site of action — the same surprising-target case we saw in the backtest.

The screen does not "know" any of these combinations explicitly — they emerge from the noisy-OR formulation of disease coverage, the requirement that the bridge target be one the primary substance does not already hit, and the target-overlap filter that rejects same-class redundancy. That known clinical combinations (BRAF + MEK in melanoma) emerge naturally as top hits in adjacent indications (BRAF + MEK in CFC syndrome) is supporting evidence that the combination finder is producing biologically meaningful candidates rather than artefacts of the scoring layer.

### 3.6 Operational outputs

For the top 5 hypotheses at the time of writing, the deal-memo generator produced 10 PDFs (5 hypotheses × 2 regions). Each is approximately 4 KB, one A4 page, with sections weaving mechanism evidence into narrative prose, IP runway commentary, market opportunity (with prevalence source attribution), KOL credentials and contact email when present, a group-pitch block describing what BioBix brings (target validation, biomarker identification, patient stratification, mechanism-of-action mining, IIT support), a proposed collaboration model (investigator-initiated trial with the substance owner; BioBix as bioinformatics analytics partner), and a signature block.

---

## 4 Discussion

### 4.1 What is novel here

The technical content of each enrichment layer in isolation is not novel: novelty subtraction via ontology hops [refs], directionality inference from genetic evidence [refs], tissue expression gating [refs], pathway-overlap analysis [refs], and PubMed-based literature priors [refs] are all individually established. The contribution of this work is integration — assembling these layers into a single pipeline with measured precision, with the regulatory and economic signals (IP, prevalence, generic status) that determine whether a candidate is actionable, with the human-network signal (identified KOLs with credentials and contact details) that determines whether a candidate is approachable, and with the operational artefact (PDF deal memos) that determines whether the work can leave the screen and reach a clinical-partnership conversation.

The backtest validation framework is also a contribution. Many repurposing screens are presented with anecdotal example hits; we instead measure precision against an explicit ground-truth set, identify the failure modes, and ship the measurement alongside the screen so any user can re-verify or extend.

### 4.2 Limitations

**Data coverage**: the validation backtest's three failures (minoxidil/alopecia, propranolol/hemangioma, acetazolamide/IIH) trace to gaps in Open Targets' disease-association data for rare or peripheral indications. The pipeline's precision is therefore bounded by OT coverage; expansion of OT (which happens with each quarterly release) raises the ceiling automatically, but novel mechanisms not yet curated in any biomedical database will not be surfaced.

**Severity heuristic narrowness**: the receptor-LoF + agonist heuristic relies on disease-name keywords ("deficiency", "complete absence", etc.) rather than structured severity annotations. A more principled solution would integrate ClinGen mechanism-of-disease assertions or OMIM clinical-feature curation, both of which are inconsistently structured at present.

**KOL disambiguation**: Semantic Scholar uses abbreviated first names that may collide on common surnames (we showed a `Boris Keren` example with five candidate authors). The highest-h-index heuristic resolves the common case correctly but cannot resolve cases where a true KOL is genuinely junior (low h-index). Affiliation-based disambiguation against the PubMed institution string would be a natural next refinement but Semantic Scholar's `affiliations` field is often empty.

**Literature pass scaling**: at the default `top_n = 5,000`, the literature pass takes ~10 minutes wall-clock with an NCBI API key. For full-scale enrichment of all 175k+ hypotheses, an ~6-hour run is required, or a separate parallel pre-cache strategy.

**IP scope**: the Orange Book covers small molecules only. Biologics IP must be sourced from the FDA Purple Book or pharma annual reports; this is a known gap that the `NA` values on biologic columns correctly signal.

**Combination-therapy scope**: the current combo finder evaluates two-drug combinations only and limits the search to companions whose target-overlap with the primary is below 50%. Three-drug combinations and more sophisticated synergy modelling (e.g., Bliss independence, Loewe additivity) are obvious extensions.

### 4.3 Comparison to existing tools

[Brief comparison to DrugCentral, OpenTargets Platform UI, PREDICT, REMEDIORE, Reactome-FI, etc.]

### 4.4 Future directions

The pipeline is engineered for incremental extension. Natural next layers:
- Cell-line / DepMap evidence for cancer hypotheses
- Structured drug-drug interaction filtering (DrugBank or FDA labels) for combination-therapy safety
- Multi-omics-based patient stratification analyses for the operational-output layer
- Integration of off-label use registries (LabRx, Medicare claims) for hypothesis grounding
- A web-deployed dashboard for non-technical end-users

A nightly refresh pipeline with version-pinned data sources (already half-engineered in the existing `download_data.sh` and `build_*` scripts) would let the screen become a continuously-updated resource rather than a periodic snapshot.

---

## 5 Code and data availability

The full pipeline source, configuration, curated CSVs, and the validation YAML are publicly available at [REPOSITORY URL TBD]. All third-party data sources are linked from `scripts/download_data.sh`. The screen runs end-to-end on a single laptop with ~50GB free disk and ~16GB RAM.

---

## 6 Acknowledgements

[To be added.]

## 7 References

[To be added. Approximately 50-80 citations covering: drug repurposing background, Open Targets, ChEMBL, Reactome, STRING, EFO, EFO+ontology novelty methods, mechanism-driven screens (PREDICT, DrugCentral), L1000/CMap, network pharmacology, individual drug repurposing successes used in the validation set, FDA Orange Book, Orphanet, PubMed, Semantic Scholar, IMPC + PhenoDigm, GBD/IHME, and methods references for noisy-OR aggregation and ontology-graph distance.]
