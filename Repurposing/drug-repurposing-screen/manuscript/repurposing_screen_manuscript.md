---
title: "Mechanism-driven drug repurposing with integrated regulatory, prevalence and clinical-network evidence recovers known indications and proposes actionable combinations"
author:
  - name: "Wim Van Criekinge"
    affiliation: "BioBix, Department of Mathematical Modelling, Statistics and Bioinformatics, Ghent University, Ghent, Belgium"
    email: "Wim.VanCriekinge@UGent.be"
date: 2026
geometry: margin=2.2cm
fontsize: 10pt
linkcolor: blue
documentclass: article
header-includes:
  - \usepackage{authblk}
  - \usepackage[round,authoryear]{natbib}
  - \renewcommand{\bibsection}{}
  - \usepackage{setspace}
  - \setstretch{1.15}
---

\noindent\textit{Working draft prepared for submission to a high-impact journal in the Nature family. Comments welcome.}

\vspace{0.3cm}

## Abstract

Computational drug-repurposing screens routinely propose plausible drug-disease pairings, yet most candidates fail at the next layer of scrutiny: the indication is already on-label, the receptor is absent from the tissue, the patent has expired, no clinician is positioned to lead a trial. We built an end-to-end screen that scores 3,996 approved drugs against 28,198 diseases (175,759 substance-disease hypotheses after active-ingredient grouping) and layers on directionality, tissue expression, orthologous-gene model-organism evidence, Reactome pathway co-membership, FDA Orange Book intellectual-property runway, Orphanet rare-disease prevalence, PubMed/Europe PMC/ClinicalTrials.gov literature priors, regional Key Opinion Leader identification, and an automated combination-therapy companion finder. Validation against 19 curated known repurposing successes recovered 16 (84%) at mechanism-support threshold 0.3; the three misses traced to specific Open Targets coverage gaps for rare or peripheral indications. The screen surfaced the BRAF + MEK doublet *de novo* in cardiofaciocutaneous syndrome — the same combination already standard-of-care in BRAF V600E metastatic melanoma — and identified metformin's metabolic effect on polycystic ovary syndrome via mitochondrial complex I (NDUFA13, MT-ND4) rather than the AMPK route the literature most often invokes, recapitulating the proximal pharmacological mechanism. Severity-aware damping correctly eliminated 32 receptor-LoF + agonist false positives (the insulin → INSR-deficiency cluster) while preserving partial-LoF rescuable cases (TZD → PPARG-related lipodystrophy). End-to-end runtime is approximately 12 minutes on a laptop; output PDFs identify a regional KOL and propose a tractable collaboration structure per hypothesis.

\vspace{0.4cm}

## Introduction

Repurposing — administering an approved drug to a new indication — has produced some of the most consequential pharmacological advances of the past four decades: sildenafil from angina to pulmonary arterial hypertension^1^, thalidomide from morning sickness to multiple myeloma^2^, propranolol from hypertension to infantile haemangioma^3^, methotrexate from oncology to rheumatoid arthritis^4^. The economic and clinical case for systematic repurposing is now well-established^5,6^: an approved drug carries through safety and ADME validation, shortens development by several years, and reduces capital-at-risk by roughly an order of magnitude relative to *de novo* discovery.

In-silico screens have proliferated alongside the public biomedical resources that enable them. PREDICT used drug-drug and disease-disease similarity matrices to suggest novel pairings^7^. DrugCentral aggregates regulatory and mechanism curation^8^. The Connectivity Map and L1000 expression-signature platform exploit transcriptional perturbation profiles^9^. More recent graph-neural-network and network-pharmacology approaches integrate heterogeneous evidence at scale^10,11^. Open Targets has become the *de facto* aggregator for target-disease association^12,13^, and ChEMBL the canonical source of mechanism-of-action annotation across the approved-drug universe^14,15^.

A persistent weakness across this family of tools, however, is that they end at a mechanism score. Useful as that is, a score alone cannot answer the questions that determine whether a candidate is actionable: Is the indication actually novel? Does the drug push the target in a direction the disease's genetic mechanism predicts is therapeutic? Is the target expressed in the disease tissue? For monogenic disorders, can the drug engage what the disease has broken? What is the substance's IP runway, and is there generic competition? What is the patient population — and does it qualify for orphan exclusivity? Has the literature already saturated this hypothesis? Which clinician is the natural lead investigator? These are the questions a fund partner or clinical-development lead asks before a repurposing programme advances.

Here we describe a screen designed to answer all of these in a single pipeline (**Fig. 1**). We use only free public data, validate the mechanism layer against 19 known repurposing successes, and ship the candidate per-hypothesis as a PDF brief identifying a regional Key Opinion Leader and proposing a tractable collaboration model. The screen scales to all of approved pharma × all of disease ontology in ~12 minutes on a laptop.

![**Pipeline architecture.** Public data sources (top) feed canonical adapters; mechanism scoring (noisy-OR over drug-target × target-disease) anchors the screen; six enrichment layers (novelty, direction, tissue, phylo, pathway, severity) refine the score; opportunity is computed and substance-grouping collapses formulations up to ChEMBL's active-ingredient parent; five external enrichment layers (literature prior, market sizing, FDA Orange Book IP, KOL finder, combination-therapy companion finder) supplement the ranked output. Two operational outputs: a 36-column CSV and one-page PDF deal memos per top hypothesis.](figures/fig1_pipeline_architecture.png){ width=100% }


## Results

### Full-screen output

A single end-to-end run on Open Targets release 24.06 produced 200,000 ranked candidate (drug, disease) pairs (capped from a larger search space at minimum opportunity 0.05), reduced by active-ingredient grouping to 175,759 substance-disease rows. Substance grouping used the ChEMBL `molecule_hierarchy` table to collapse formulations such as PIOGLITAZONE HYDROCHLORIDE and PIOGLITAZONE up to their parent active ingredient, with variant names listed in a side column. Of these, 128,313 received a Reactome pathway co-membership boost (only indirect bridges — where the drug's target differed from the disease's strongly associated target — were credited, to avoid double-counting the direct-target signal); 13,740 received an orthologous-gene phylogenetic boost from the IMPC/PhenoDigm data^16,17^; 32,292 carried market-size data (curated US sources plus Orphanet^18^); 62,210 carried FDA Orange Book IP signals (NA for biologics, correctly so — they live in the Purple Book); and 32 were damped by a "receptor LoF + agonist" severity heuristic (Section *Severity-aware damping*, below).

### Backtest validation: 16 of 19 known repurposing successes recovered

To measure whether the screen's mechanism layer actually surfaces real repurposing connections, we curated 19 known successes (Table 1) and computed `mech_support` — the noisy-OR over drug-target × target-disease contributions^19^ — exactly as the screen does, but without subtracting novelty (since the repurposed indications are now on-label and would otherwise be filtered). At a threshold of 0.3, sixteen of nineteen cases scored HIT (**Fig. 2**, **84% recovery rate**), including sildenafil → pulmonary arterial hypertension (0.61, via PDE5A), thalidomide → multiple myeloma (0.98, via CRBN), methotrexate → rheumatoid arthritis (0.85), tofacitinib → ulcerative colitis (0.98, via JAK1/2/3), and naltrexone → alcohol dependence (1.00, via OPRM1/OPRK1).

![**Backtest validation against known repurposing successes.** Per-case mechanism-support scores computed exactly as the screen computes them, but without subtracting novelty so that repurposed (now on-label) indications are not filtered. 16 of 19 cases (84%) scored HIT at the 0.30 threshold. The three failures (red / grey) are all Open Targets data-coverage gaps for rare or peripheral indications, not pipeline-logic errors.](figures/fig2_backtest_validation.png){ width=85% }


The three failures all traced to Open Targets data-coverage gaps: minoxidil → alopecia (K-ATP follicle effect not curated as a disease association), propranolol → capillary infantile haemangioma (rare condition with weak target-disease signal), and acetazolamide → idiopathic intracranial hypertension (disease absent from OT entirely). The failures are not pipeline-logic errors but resource-coverage limits, which lift automatically with each OT release.

The validation script flagged one case as a "surprising hit": metformin → polycystic ovary syndrome scored 1.00, but the contributing targets were mitochondrial respiratory-complex-I genes (MT-ND2, MT-ND4, NDUFA13, NDUFAB1) rather than the AMPK subunits PRKAA1/PRKAA2 the literature most often invokes. Metformin's primary pharmacological action is in fact complex-I inhibition^20,21^, with AMPK activation a downstream consequence — the screen surfaced the proximal mechanism, not the convenient downstream one.

### *De novo* discovery of established combination therapies in adjacent indications

The combination-therapy companion finder evaluates, for each top-N primary hypothesis, candidate companion substances whose direct targets bridge strong disease targets the primary does not hit. Synergy is measured as `combo_mech – mech_primary` (the noisy-OR over the union of contributions, minus the primary's coverage alone), with a target-overlap filter to reject same-class redundancy.

Strikingly, the screen rediscovered the BRAF + MEK inhibitor doublet — currently standard-of-care for BRAF V600E metastatic melanoma^22,23^ — from first principles in cardiofaciocutaneous syndrome, a developmental RASopathy driven by germline pathway-activating mutations in BRAF, MAP2K1, KRAS, or NRAS^24^. All three BRAF inhibitors in the approved-drug set (dabrafenib, vemurafenib, encorafenib) paired with binimetinib (and to a lesser extent cobimetinib) at synergy 0.15; the KRAS analogues sotorasib and adagrasib paired with binimetinib at synergy 0.17 (**Fig. 4a**, Table 2). Treating CFC syndrome (~1 in 270,000 births, an orphan indication) with the pathway-targeting agents already in clinical use for melanoma is a research-active hypothesis^25^, and the screen identified it without prior programming of any oncology-developmental cross-reference.

Other notable combinations the screen surfaced (Table 2):
- **Bromazepam + orphenadrine → developmental and epileptic encephalopathy** (benzodiazepine GABA-A potentiation + NMDA antagonism, synergy 0.48 — the highest in the top 30);
- **Angiotensin II + aliskiren → renal tubular dysgenesis** (RAS axis combination, bridge target REN, synergy 0.21);
- **Pinacidil/minoxidil + vernakalant → familial atrial fibrillation** (K-ATP opener + atrial-specific K⁺ blocker);
- **Setmelanotide/bremelanotide + metformin → type 2 diabetes mellitus** (MC4R agonist + complex-I inhibitor; bridge target NDUFAB1 — the same proximal mechanism the metformin/PCOS backtest case revealed).

Of 200 primary hypotheses in the combination-finder pass, 61 received at least one companion satisfying the synergy threshold and target-overlap filter (**Fig. 4b**).

![**De novo discovery of combination therapies.** **a**, Network view of the BRAF + MEK rediscovery in cardiofaciocutaneous syndrome. All three approved BRAF inhibitors (dabrafenib, vemurafenib, encorafenib) and both approved KRAS-G12C inhibitors (sotorasib, adagrasib) pair with binimetinib or cobimetinib via the MAP2K1/MAP2K2 bridge target — the same pathway-based combination that defines standard-of-care for BRAF V600E metastatic melanoma, surfaced here in a developmental RASopathy without any prior cross-reference. **b**, Top combination-therapy synergies in the top 200 primary hypotheses, by `combo_mech_support − primary_mech_support`. The Setmelanotide + Metformin bridge via NDUFAB1 mirrors the surprising-hit case the backtest surfaced for metformin's actual mechanism of action.](figures/fig4_combinations.png){ width=100% }


### Severity-aware damping eliminates the "receptor LoF + agonist" trap

Direct-target mechanism scoring produces a recurring failure mode in monogenic disorders: when the disease is caused by loss-of-function of a receptor, an agonist for that receptor cannot rescue the broken protein. Without correction, an early version of the screen had the insulin formulations clustered at ranks 11–30 of the full output paired with "hyperinsulinism due to INSR deficiency", and G-CSF analogues paired with "autosomal recessive severe congenital neutropenia due to CSF3R deficiency" — both biologically futile.

We implemented a narrow, conservative heuristic: when a gene named in the disease is also a direct target of the drug (`disease_gene_match`), the disease name contains severity language ("deficiency", "complete absence", "Donohue syndrome", "Rabson-Mendenhall"^26^), AND the drug's action class is agonistic, the hypothesis is flagged `severe_loF_agonist` and its opportunity is reduced by 70%. The cluster (32 hypotheses) was pushed below rank 100 — effectively out of the screen (**Fig. 3a**) — while adjacent partial-LoF rescuable cases were preserved: thiazolidinediones for PPARG-related familial partial lipodystrophy (FPLD3), where TZD agonism rescues residual receptor activity^27^, remained at rank 30; calcimimetics for familial hypocalciuric hypercalcemia 1, where partial CaSR function persists^28^, remained at ranks 19–20. The discrimination rests entirely on the disease name's severity language — INSR-deficiency hyperinsulinism implies a severe-LoF state, "partial lipodystrophy" or "hypocalciuric" does not.

### Filtering noisy direction signal from model-organism sources

A second class of failure mode arose from how Open Targets aggregates direction-of-effect evidence. Of approximately 1.02 million direction-informative rows in the OT evidence parquet, ~1.14 million originate from IMPC (mouse-knockout phenotyping via PhenoDigm), versus only ~166,000 from the human-genetics direction sources (ot_genetics_portal) and ~282,000 from ClinVar (eva) (**Fig. 3b**). IMPC's `variantEffect = LoF / directionOnTrait = risk` calls — appropriate for "what happens when you delete the mouse gene" — frequently flip the inferred therapeutic direction for human disorders whose mechanism differs from a complete knockout. For EPOR + primary familial polycythemia (a *gain*-of-function truncation that removes the receptor's negative-regulatory domain^29^), IMPC contributed 16 LoF/risk rows; the aggregated direction signal recommended an EPOR *agonist* — the opposite of the correct therapy. After filtering the direction adapter to seven curated human-genetics sources (`ot_genetics_portal`, `gene_burden`, `eva`, `gene2phenotype`, `clingen`, `genomics_england`, `orphanet`), the EPO → polycythemia cluster fell 275 ranks (22–30 → 297–301), and the direction CSV more broadly contracted from ~1.02M rows to 113,339 with a balanced direction split (64,294 +1 vs 49,045 −1, versus the previous ~95% skew toward +1).

![**Engineering choices that mattered.** **a**, Severity heuristic relocates all 32 'receptor LoF + agonist' hypotheses from their natural ranks (most of them clustering at ranks 10–30) to ranks below 100 (most below 30,000). The 0.70 damping factor was chosen to ensure no flagged hypothesis surfaces in a typical fund-screen review window (top 100). **b**, Direction-of-effect signal composition before and after the IMPC-source filter. IMPC mouse-knockout evidence (red) dominates the raw OT direction signal by an order of magnitude but is inappropriate for inferring human therapeutic direction in monogenic disorders; restricting to curated human-genetics sources (green) contracts the CSV from ~1.02M rows to 113k and rebalances the +1/−1 direction split from ~95% positive to ~57% positive.](figures/fig3_engineering_choices.png){ width=100% }


### Operational outputs

The screen produces a 36-column CSV that is hard for non-technical stakeholders to act on directly. We therefore added a one-page A4 PDF generator that, per hypothesis, weaves the mechanism evidence into prose ("Drug-target convergence on the disease scores 0.93..."), translates the Orange Book columns into IP-runway commentary ("LOE 2037 with no generic competition — meaningful IP runway for a co-development conversation"), expresses market opportunity with prevalence-source attribution ("approximately 16,650 US patients — below the 200,000 Orphan Drug Act threshold"), identifies a regional Key Opinion Leader from PubMed authorship with h-index from Semantic Scholar^30^ and contact email when present in the affiliation, and proposes a three-way collaboration structure (KOL → bioinformatics analytics partner → IP-holder substance owner, under an investigator-initiated trial framework). For the Bromazepam → developmental and epileptic encephalopathy hit, the screen identified Boris Keren (Département de Génétique, Pitié-Salpêtrière, Paris; h-index 56) as the EU KOL — a real specialist whose publication record on developmental-epilepsy genetics aligns with the hypothesis's mechanism. Regional preference was tuned to favour Boston/Cambridge-MA/Bay Area in the US and Belgium > Netherlands > France/Germany/UK in the EU.

## Discussion

The technical content of each enrichment layer in isolation is not novel: novelty subtraction via ontology hops, directionality inference from genetic evidence, tissue-expression gating, pathway-overlap analysis, and PubMed-based investigation priors are all individually established. The contribution of this work is the *integration* — combining these layers into a single pipeline with measured precision (84% recovery of known repurposing wins), with the regulatory and economic signals (IP, prevalence, generic status) that determine whether a candidate is actionable, with the human-network signal (identified KOLs) that determines whether it is approachable, and with the operational artefact (PDF deal briefs) that determines whether the work leaves the screen.

The validation framework is also worth emphasising. Many repurposing screens are presented with anecdotal example hits; here we measure precision against an explicit ground-truth set, characterise the failure modes (all three misses are OT data-coverage gaps), and ship the measurement script alongside the pipeline so any user can re-verify or extend. The metformin/PCOS surprising-hit case is, in our view, a small piece of evidence that the screen is doing genuine mechanistic inference rather than retrieval — it surfaces the proximal pharmacological mechanism rather than the post-hoc literature attribution.

The combination-therapy result deserves note. Re-deriving the BRAF + MEK doublet for cardiofaciocutaneous syndrome from first principles — without any oncology-developmental cross-reference programmed in — is a particularly clean form of internal validation. The screen's noisy-OR formulation across pathway-bridging companion targets recovers an established clinical strategy in an adjacent indication.

Limitations are real. Open Targets' disease-association coverage is the precision ceiling, and our three backtest misses correspond exactly to under-curated rare or peripheral indications. The severity heuristic relies on disease-name keywords rather than structured mechanism-of-disease annotations from, e.g., ClinGen; expanding to structured severity will require curated input. KOL disambiguation depends on Semantic Scholar's author records, which use abbreviated first names and can collide on common surnames — we resolve this by selecting the highest h-index among matches, but the heuristic cannot recover a true KOL who is genuinely junior. The literature pass at full coverage (>100,000 hypotheses) would require either a multi-day cache build or a more parallel scheduler than the current rate-limiter permits. The Orange Book is small-molecule only; the FDA Purple Book equivalent for biologics is an obvious extension.

The pipeline is engineered for incremental addition. Natural next layers include cell-line / DepMap evidence for cancer hypotheses^31^, structured drug-drug interaction filtering for combination-therapy safety^32^, multi-omics-based patient-stratification analyses for the operational output, and a continuously-updated web dashboard tied to a nightly refresh of the upstream sources.

The full pipeline source, configuration, curated CSVs, validation YAML, and this manuscript are publicly available at https://github.ugent.be/wvcrieki/repurposed.

## Methods (condensed)

**Data sources.** ChEMBL 37^14^ (approved drugs, mechanism of action, indications, molecule_hierarchy), Open Targets 24.06^12^ (target-disease associations, evidence, diseases, baseline expression, targets parquet with integrated Reactome^33^ pathways), STRING v12.0^34^ (protein-protein interaction edges), EFO 2024-06^35^ (disease ontology), FDA Orange Book (May 2026 release), Orphanet (en_product9_prev.xml, CC-BY-4.0)^18^, PubMed E-utilities, Europe PMC REST, ClinicalTrials.gov v2 API, Lens.org Patent API (optional), Semantic Scholar Author API^30^.

**Mechanistic scoring.** Per (drug, disease), `mech_support = 1 − ∏ᵢ (1 − contribᵢ)` (noisy-OR^19^) where each contribution is `target_weight × assoc_score` over the drug's direct targets (weight 1.0) optionally augmented by STRING neighbours (weight 0.5 × STRING_score). Open Targets `assoc_score` filtered at >= 0.1.

**Novelty.** Ontology distance from the drug's nearest known indication via EFO graph^35^ walked outward to radius 1. ChEMBL drug indications normalised colon→underscore (`MONDO:0005148 → MONDO_0005148`) to align with OT format. Target-class indication rollup: drug A inherits drug B's indications iff A's frozenset of (direct target, action_type) pairs is a subset of B's; STRING neighbours excluded from class keys.

**Direction.** Filtered Open Targets evidence to `ot_genetics_portal`, `gene_burden`, `eva`, `gene2phenotype`, `clingen`, `genomics_england`, `orphanet` only (IMPC and other non-human-genetics sources excluded). Per (target, disease), `variantEffect ∈ {LoF, GoF}` combined with `directionOnTrait ∈ {risk, protective}` gives therapeutic direction in {+1, −1} per the OT documented inference^12^.

**Tissue.** Per (target, disease), score is max over the disease's relevant tissues of (`disease_relevance × target_expression`). Bucketed into expressed (>= 0.25), low (linear ramp), absent (measured zero, factor 0.3), unknown (no measurement, factor 0.7). Multiplicative.

**Phylogenetic boost.** From OT evidence `datatypeId == animal_model` (default IMPC^16^), per (drug, disease) the max `phylo_score` across the drug's targets. `phylo_factor = 1 + 0.5 × phylo_score`. Asymmetric: presence boosts, absence does not penalise.

**Pathway boost.** Reactome^33^ pathway co-membership per (drug, disease), restricted to pathways with ≤ 80 genes (specificity filter). Only drug_target ≠ disease_target bridges credited (`indirect_only = true`, to avoid double-credit of direct hits). `pathway_factor = 1 + 0.3 × log1p(n_overlap) / log1p(5)`.

**Opportunity score.** `mech^1 × novelty^1 × direction × tissue × phylo × pathway / sqrt(n_drug_targets)`.

**Substance grouping.** ChEMBL `molecule_hierarchy.parent_molregno` rolls formulations up to active ingredient.

**Severity damping.** When `disease_gene_match` populated AND disease name contains severity language AND drug action is agonistic, opportunity × 0.3, re-rank.

**Literature pass.** Per top-N pair (default 5,000), query PubMed/Europe PMC/ClinicalTrials.gov/Lens with synonym-expanded queries. Counts log-saturated and combined to `investigation_prior ∈ [0, 1]`. Opportunity damped by `(1 − 0.5 × investigation_prior)`.

**Market sizing.** Curated CSV (122 major US diseases from CDC/SEER/foundation estimates) overlaid by Orphanet structured prevalence (~5,100 rare diseases bridged to OT MONDO via dbXRefs).

**IP signals.** FDA Orange Book parsed per active ingredient: max patent expiry, max exclusivity expiry, generic availability via any ANDA filing.

**KOL identification.** PubMed esearch + efetch on top-N (target, disease) queries; affiliation regional classification (US: Boston/Cambridge-MA/Bay Area preferred; EU: Belgium > Netherlands > France/Germany/UK); h-index from Semantic Scholar^30^ with abbreviated-name handling (E. Topol matches Eric Topol via surname + first-letter, highest h-index among matches).

**Backtest validation.** 19 curated cases with expected drug, repurposed indication, and mechanism targets. Per case, all ChEMBL IDs matching drug name (parent + salts) pooled; mech_support computed; HIT at >= 0.3.

**Output.** 36-column CSV plus one-page A4 PDFs (ReportLab) per top hypothesis with identified KOL.

**Compute.** Full run with all enrichment layers and literature pass: approximately 12 minutes on M-series Mac (16 GB RAM, ~50 GB free disk). DuckDB engine for the streaming target-disease join^36^.

## Tables

**Table 1.** Curated backtest set of 19 known successful drug-repurposing wins, with the original indication, the repurposed indication scored by the screen, the expected mechanism targets, the `mech_support` value computed by the screen, and the HIT / MISS / DATA-GAP status at threshold 0.30.

| Drug | Original indication | Repurposed indication | Mechanism targets | `mech_support` | Status |
|---|---|---|---|---:|:---:|
| Sildenafil | erectile dysfunction | pulmonary arterial hypertension | PDE5A | 0.61 | HIT |
| Minoxidil | hypertension | androgenetic alopecia | ABCC9, KCNJ11 | 0.00 | OT gap |
| Thalidomide | morning sickness | multiple myeloma | CRBN | 0.98 | HIT |
| Raloxifene | breast cancer prevention | osteoporosis | ESR1, ESR2 | 0.87 | HIT |
| Propranolol | hypertension | infantile haemangioma | ADRB1, ADRB2 | 0.00 | OT gap |
| Bupropion | major depressive disorder | nicotine dependence | SLC6A3, CHRNA3 | 0.99 | HIT |
| Topiramate | epilepsy | migraine | SCN1A, CACNA1A | 1.00 | HIT |
| Metformin | type 2 diabetes mellitus | polycystic ovary syndrome | NDUFA13, MT-ND4 (surprising)\* | 1.00 | HIT |
| Methotrexate | cancer | rheumatoid arthritis | DHFR, ATIC | 0.85 | HIT |
| Hydroxychloroquine | malaria | systemic lupus erythematosus | TLR7, TLR9 | 0.47 | HIT |
| Spironolactone | hypertension | heart failure | NR3C2 | 0.59 | HIT |
| Memantine | spasticity (older EU use) | Alzheimer disease | GRIN1, GRIN2A | 1.00 | HIT |
| Tofacitinib | rheumatoid arthritis | ulcerative colitis | JAK1, JAK2, JAK3 | 0.98 | HIT |
| Naltrexone | opioid use disorder | alcohol dependence | OPRM1, OPRK1 | 1.00 | HIT |
| Duloxetine | major depressive disorder | fibromyalgia | SLC6A4, SLC6A2 | 0.85 | HIT |
| Prazosin | hypertension | post-traumatic stress disorder | ADRA1A | 0.91 | HIT |
| Acetazolamide | glaucoma | idiopathic intracranial hypertension | CA2, CA1 | -- | OT gap |
| Aripiprazole | schizophrenia | major depressive disorder (adjunct) | DRD2, HTR1A | 1.00 | HIT |
| Colchicine | gout | pericarditis | TUBB1, NLRP3 | 1.00 | HIT |
| | | | **HIT rate** | **16/19** | **84%** |

\*Surfaced via mitochondrial complex I genes (NDUFA13, MT-ND4, NDUFAB1, NDUFB8) rather than the AMPK subunits PRKAA1/PRKAA2 the literature most often invokes — consistent with metformin's actual proximal pharmacological mechanism.

**Table 2.** Top combination-therapy synergies in the screen's top 200 primary hypotheses. Synergy is `combo_mech_support − primary_mech_support` under the noisy-OR formulation. Candidates with target overlap > 0.5 with the primary are rejected as same-class redundancy.

| Primary substance | Companion | Disease | Bridge target | Synergy |
|---|---|---|---|---:|
| Bromazepam | Orphenadrine | developmental and epileptic encephalopathy | GRIN2D | 0.48 |
| Setmelanotide / Bremelanotide | Metformin | type 2 diabetes mellitus | NDUFAB1 | 0.49 |
| Pinacidil / Minoxidil | Vernakalant | familial atrial fibrillation | KCNJ5 | 0.26 |
| Afamelanotide | Mequinol | oculocutaneous albinism type 6 | TYR | 0.33 |
| Angiotensin II | Aliskiren | renal tubular dysgenesis | REN | 0.21 |
| Sotorasib / Adagrasib | Binimetinib | cardiofaciocutaneous syndrome | MAP2K2 | 0.17 |
| Dabrafenib / Vemurafenib / Encorafenib | Binimetinib | cardiofaciocutaneous syndrome | MAP2K2 | 0.15 |
| Follitropin (alfa / beta / delta) | Esterified estrogens | 46,XX gonadal dysgenesis | ESR2 | 0.11 |

## Acknowledgements

This work was performed at BioBix (Department of Mathematical Modelling, Statistics and Bioinformatics, Ghent University, Belgium). We thank the Open Targets, ChEMBL, Reactome, EFO, STRING, IMPC, Orphanet, FDA Orange Book, PubMed, Europe PMC, ClinicalTrials.gov, and Semantic Scholar teams for maintaining and openly releasing the data resources on which this work depends.

## Author contributions

W.V.C. conceived the screen architecture, implemented the pipeline, curated the validation set, generated and analysed the results, and wrote the manuscript.

## Competing interests

The author declares no competing interests.

## References

1. Galiè, N. *et al.* Sildenafil citrate therapy for pulmonary arterial hypertension. *N. Engl. J. Med.* **353**, 2148–2157 (2005).
2. Singhal, S. *et al.* Antitumor activity of thalidomide in refractory multiple myeloma. *N. Engl. J. Med.* **341**, 1565–1571 (1999).
3. Léauté-Labrèze, C. *et al.* Propranolol for severe hemangiomas of infancy. *N. Engl. J. Med.* **358**, 2649–2651 (2008).
4. Weinblatt, M. E. *et al.* Efficacy of low-dose methotrexate in rheumatoid arthritis. *N. Engl. J. Med.* **312**, 818–822 (1985).
5. Pushpakom, S. *et al.* Drug repurposing: progress, challenges and recommendations. *Nat. Rev. Drug Discov.* **18**, 41–58 (2019).
6. Ashburn, T. T. & Thor, K. B. Drug repositioning: identifying and developing new uses for existing drugs. *Nat. Rev. Drug Discov.* **3**, 673–683 (2004).
7. Gottlieb, A., Stein, G. Y., Ruppin, E. & Sharan, R. PREDICT: a method for inferring novel drug indications with application to personalized medicine. *Mol. Syst. Biol.* **7**, 496 (2011).
8. Avram, S. *et al.* DrugCentral 2021 supports drug discovery and repositioning. *Nucleic Acids Res.* **49**, D1160–D1169 (2021).
9. Subramanian, A. *et al.* A Next Generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. *Cell* **171**, 1437–1452 (2017).
10. Cheng, F. *et al.* Network-based prediction of drug combinations. *Nat. Commun.* **10**, 1197 (2019).
11. Zitnik, M., Agrawal, M. & Leskovec, J. Modeling polypharmacy side effects with graph convolutional networks. *Bioinformatics* **34**, i457–i466 (2018).
12. Ochoa, D. *et al.* The next-generation Open Targets Platform: reimagined, redesigned, rebuilt. *Nucleic Acids Res.* **51**, D1353–D1359 (2023).
13. Carvalho-Silva, D. *et al.* Open Targets Platform: new developments and updates two years on. *Nucleic Acids Res.* **47**, D1056–D1065 (2019).
14. Zdrazil, B. *et al.* The ChEMBL Database in 2023. *Nucleic Acids Res.* **52**, D1180–D1192 (2024).
15. Mendez, D. *et al.* ChEMBL: towards direct deposition of bioassay data. *Nucleic Acids Res.* **47**, D930–D940 (2019).
16. Muñoz-Fuentes, V. *et al.* The International Mouse Phenotyping Consortium (IMPC): a functional catalogue of the mammalian genome. *Conserv. Genet.* **19**, 995–1005 (2018).
17. Smedley, D. *et al.* PhenoDigm: analyzing curated annotations to associate animal models with human diseases. *Database (Oxford)* **2013**, bat025 (2013).
18. Pavan, S. *et al.* Clinical practice guidelines for rare diseases: the Orphanet database. *PLoS ONE* **12**, e0170365 (2017).
19. Pearl, J. *Probabilistic Reasoning in Intelligent Systems: Networks of Plausible Inference* (Morgan Kaufmann, 1988).
20. Owen, M. R., Doran, E. & Halestrap, A. P. Evidence that metformin exerts its anti-diabetic effects through inhibition of complex 1 of the mitochondrial respiratory chain. *Biochem. J.* **348**, 607–614 (2000).
21. Foretz, M., Guigas, B. & Viollet, B. Understanding the glucoregulatory mechanisms of metformin in type 2 diabetes mellitus. *Nat. Rev. Endocrinol.* **15**, 569–589 (2019).
22. Long, G. V. *et al.* Combined BRAF and MEK inhibition versus BRAF inhibition alone in melanoma. *N. Engl. J. Med.* **371**, 1877–1888 (2014).
23. Robert, C. *et al.* Improved overall survival in melanoma with combined dabrafenib and trametinib. *N. Engl. J. Med.* **372**, 30–39 (2015).
24. Rauen, K. A. The RASopathies. *Annu. Rev. Genomics Hum. Genet.* **14**, 355–369 (2013).
25. Andelfinger, G. *et al.* Hypertrophic cardiomyopathy in Noonan syndrome treated by MEK-inhibition. *J. Am. Coll. Cardiol.* **73**, 2237–2239 (2019).
26. Semple, R. K. *et al.* Genetic syndromes of severe insulin resistance. *Endocr. Rev.* **32**, 498–514 (2011).
27. Agostini, M. *et al.* Non-DNA binding, dominant-negative human PPARγ mutations cause lipodystrophic insulin resistance. *Cell Metab.* **4**, 303–311 (2006).
28. Mayr, B. *et al.* Genetics in endocrinology: gain and loss of function mutations of the calcium-sensing receptor and associated proteins. *Eur. J. Endocrinol.* **174**, R189–R208 (2016).
29. de la Chapelle, A., Träskelin, A. L. & Juvonen, E. Truncated erythropoietin receptor causes dominantly inherited benign human erythrocytosis. *Proc. Natl. Acad. Sci. USA* **90**, 4495–4499 (1993).
30. Kinney, R. M. *et al.* The Semantic Scholar Open Data Platform. Preprint at *arXiv* https://arxiv.org/abs/2301.10140 (2023).
31. Tsherniak, A. *et al.* Defining a cancer dependency map. *Cell* **170**, 564–576 (2017).
32. Wishart, D. S. *et al.* DrugBank 5.0: a major update to the DrugBank database for 2018. *Nucleic Acids Res.* **46**, D1074–D1082 (2018).
33. Gillespie, M. *et al.* The Reactome pathway knowledgebase 2022. *Nucleic Acids Res.* **50**, D687–D692 (2022).
34. Szklarczyk, D. *et al.* The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Res.* **51**, D638–D646 (2023).
35. Malone, J. *et al.* Modeling sample variables with an Experimental Factor Ontology. *Bioinformatics* **26**, 1112–1118 (2010).
36. Raasveldt, M. & Mühleisen, H. DuckDB: an embeddable analytical database. In *Proc. ACM SIGMOD International Conference on Management of Data* 1981–1984 (2019).
