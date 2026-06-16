---
title: "REPRISE: a Repurposing Engine for Pathway-Resolved Indication Scoring and Evidence — recovers known indications and proposes actionable combinations"
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

Computational drug-repurposing screens routinely propose plausible drug-disease pairings, yet most candidates fail at the next layer of scrutiny: the indication is already on-label, the receptor is absent from the tissue, the patent has expired, the underlying biology is dose-impractical. We built an end-to-end screen that scores 3,996 approved drugs against 28,198 diseases (176,272 substance-disease hypotheses after active-ingredient grouping) and layers on directionality, tissue expression, orthologous-gene model-organism evidence, Reactome pathway co-membership, FDA Orange Book intellectual-property runway, Orphanet rare-disease prevalence, PubMed/Europe PMC/ClinicalTrials.gov/Lens.org literature and patent priors, and an automated combination-therapy companion finder. Validation against 54 curated known repurposing successes recovered 48 (89%) at mechanism-support threshold 0.3; the six misses traced to specific Open Targets coverage gaps for rare or peripheral indications. The screen surfaced the BRAF + MEK doublet *de novo* in cardiofaciocutaneous syndrome — the same combination already standard-of-care in BRAF V600E metastatic melanoma — and identified metformin's metabolic effect on polycystic ovary syndrome via mitochondrial complex I (NDUFA13, MT-ND4) rather than the AMPK route the literature most often invokes, recapitulating the proximal pharmacological mechanism. Severity-aware damping correctly eliminated 32 receptor-LoF + agonist false positives (the insulin → INSR-deficiency cluster) while preserving partial-LoF rescuable cases (TZD → PPARG-related lipodystrophy). End-to-end runtime is approximately 12 minutes on a laptop; each top hypothesis ships as a one-page PDF brief that weaves the mechanism, regulatory and market evidence into a narrative suitable for clinical-collaborator outreach.

\vspace{0.4cm}

## Introduction

Repurposing — administering an approved drug to a new indication — has produced some of the most consequential pharmacological advances of the past four decades: sildenafil from angina to pulmonary arterial hypertension^1^, thalidomide from morning sickness to multiple myeloma^2^, propranolol from hypertension to infantile haemangioma^3^, methotrexate from oncology to rheumatoid arthritis^4^. Across a curated set of 54 such landmark wins assembled for this work (FDA / EMA approvals and well-established off-label uses, Methods), the median lag from original-use approval to repurposed-use approval is 8 years (range 1–44, excluding pre-modern drugs), and just five source indications — hypertension (7 wins), epilepsy (4), type 2 diabetes (4), rheumatoid arthritis (4), and schizophrenia (3) — supplied over 40% of all repurposings. The economic and clinical case for systematic repurposing is now well-established^5,6^: an approved drug carries through safety and ADME validation, shortens development by several years, and reduces capital-at-risk by roughly an order of magnitude relative to *de novo* discovery. Yet the historic process has been almost entirely clinician-driven, hub-indication-bound, and serendipitous — propranolol → haemangioma took 44 years from its hypertension approval; methotrexate → rheumatoid arthritis took 35; metformin → polycystic ovary syndrome remains off-label despite three decades of mechanistic and clinical evidence. The combinatorial search space — every approved drug against every disease — is 176,272-large in a single Open Targets release; only a sliver of it has been examined.

In-silico screens have proliferated alongside the public biomedical resources that enable them. PREDICT used drug-drug and disease-disease similarity matrices to suggest novel pairings^7^. DrugCentral aggregates regulatory and mechanism curation^8^. The Connectivity Map and L1000 expression-signature platform exploit transcriptional perturbation profiles^9^. More recent graph-neural-network and network-pharmacology approaches integrate heterogeneous evidence at scale^10,11^. Open Targets has become the *de facto* aggregator for target-disease association^12,13^, and ChEMBL the canonical source of mechanism-of-action annotation across the approved-drug universe^14,15^.

A persistent weakness across this family of tools, however, is that they end at a mechanism score. Useful as that is, a score alone cannot answer the questions that determine whether a candidate is actionable: Is the indication actually novel? Does the drug push the target in a direction the disease's genetic mechanism predicts is therapeutic? Is the target expressed in the disease tissue? For monogenic disorders, can the drug engage what the disease has broken? What is the substance's IP runway, and is there generic competition? What is the patient population — and does it qualify for orphan exclusivity? Has the literature already saturated this hypothesis, or is it commercially open? These are the questions a clinical-development lead or fund partner asks before a repurposing programme advances.

Here we describe REPRISE — a *Repurposing Engine for Pathway-Resolved Indication Scoring and Evidence* — designed to answer all of these in a single pipeline (**Fig. 1**). We use only free public data, validate the mechanism layer against 54 known repurposing successes, and ship each top candidate as a one-page PDF brief that translates the screen's numeric output into a narrative suitable for clinical-collaborator outreach. The screen scales to all of approved pharma × all of disease ontology in ~12 minutes on a laptop — compressing what has historically been a multi-year, hub-indication-bound serendipity loop into a recomputable pass with measured precision.

![**Pipeline architecture.** Public data sources (top) feed canonical adapters; mechanism scoring (noisy-OR over drug-target × target-disease) anchors the screen; six enrichment layers (novelty, direction, tissue, phylo, pathway, severity) refine the score; opportunity is computed and substance-grouping collapses formulations up to ChEMBL's active-ingredient parent; four external enrichment layers (literature and patent prior, market sizing from curated and Orphanet sources, FDA Orange Book intellectual-property runway, combination-therapy companion finder) supplement the ranked output. Two operational outputs: a 36-column CSV and per-hypothesis PDF briefs that translate the mechanism, regulatory and market evidence into a narrative for clinical-collaborator outreach.](figures/fig1_pipeline_architecture.png){ width=100% }


## Results

### Full-screen output

A single end-to-end run on Open Targets release 24.06 produced 176,272 ranked substance-disease hypotheses (active-ingredient grouped, minimum opportunity 0.05). Substance grouping used the ChEMBL `molecule_hierarchy` table to collapse formulations such as PIOGLITAZONE HYDROCHLORIDE and PIOGLITAZONE up to their parent active ingredient, with variant names listed in a side column. Of these, 128,313 received a Reactome pathway co-membership boost (only indirect bridges — where the drug's target differed from the disease's strongly associated target — were credited, to avoid double-counting the direct-target signal); 13,740 received an orthologous-gene phylogenetic boost from the IMPC/PhenoDigm data^16,17^; 32,292 carried market-size data (curated US sources plus Orphanet^18^); 62,210 carried FDA Orange Book IP signals (NA for biologics, correctly so — they live in the Purple Book); and 32 were damped by a "receptor LoF + agonist" severity heuristic (Section *Severity-aware damping*, below).

### Historical context: the pace REPRISE compresses

Across the curated 54-case set, the median lag from original-use approval to repurposed-use approval is **8 years**, but the distribution is heavy-tailed: 17 of 54 (31%) were repurposed within five years (almost all target-driven and post-2000 — imatinib → gastrointestinal stromal tumor in 1 year via KIT, semaglutide → obesity in 4 years via GLP1R, eculizumab → atypical haemolytic uremic syndrome in 4 years via C5), while 16 of 54 (30%) took 15 years or more (propranolol → infantile haemangioma 44 years; methotrexate → rheumatoid arthritis 35; thalidomide → multiple myeloma 41, after a 30-year regulatory withdrawal). The source-indication distribution is similarly concentrated: hypertension alone supplied 7 of 54 — propranolol → haemangioma and PTSD; eplerenone, spironolactone and carvedilol → heart failure; verapamil → migraine; minoxidil → alopecia — with epilepsy (4), type 2 diabetes (4), rheumatoid arthritis (4), and schizophrenia (3) together raising the share to over 40%. REPRISE re-evaluates the full equivalent screening pass — every approved drug against every disease in the public ontology, 176,272 candidate pairs — in approximately 12 minutes, with mechanism, directionality, tissue, pathway, IP-runway, prevalence, and literature-density overlays applied throughout. The public dashboard companion (Supplementary Fig. 2) visualises this historical pattern as a disease-disease network where hub indications are nodes with high degree, alongside the screen's full ranked output.

### Backtest validation: 48 of 54 known repurposing successes recovered

To measure whether the screen's mechanism layer actually surfaces real repurposing connections, we curated 54 known successes (Table 1; FDA / EMA approvals and well-established off-label uses across oncology, cardiometabolic, neuro / psychiatry, immunology, and rare disease) and computed `mech_support` — the noisy-OR over drug-target × target-disease contributions^19^ — exactly as the screen does, but without subtracting novelty (since the repurposed indications are now on-label and would otherwise be filtered). At a threshold of 0.3, forty-eight of fifty-four cases scored HIT (**Fig. 2**, **89% recovery rate**, 95% Clopper–Pearson CI 77–96%), including sildenafil → pulmonary arterial hypertension (0.61, via PDE5A), thalidomide → multiple myeloma (0.98, via CRBN), methotrexate → rheumatoid arthritis (0.85), tofacitinib → ulcerative colitis (0.98, via JAK1/2/3), naltrexone → alcohol dependence (1.00, via OPRM1/OPRK1), imatinib → gastrointestinal stromal tumor (via KIT/PDGFRA), empagliflozin → heart failure (via SLC5A2), semaglutide → obesity (via GLP1R), and lamotrigine → bipolar disorder (1.00, via SCN1A/2A).

![**Backtest validation against known repurposing successes.** Per-case mechanism-support scores computed exactly as the screen computes them, but without subtracting novelty so that repurposed (now on-label) indications are not filtered. 48 of 54 cases (89%) scored HIT at the 0.30 threshold (Clopper–Pearson 95% CI 77–96%). The six failures (red / grey) are Open Targets data-coverage gaps for rare or peripheral indications, not pipeline-logic errors.](figures/fig2_backtest_validation.png){ width=85% }


The six failures all traced to Open Targets data-coverage gaps for indications where the canonical mechanism target is not in OT's target-disease curation: minoxidil → alopecia (K-ATP follicle effect not curated), propranolol → capillary infantile haemangioma (rare condition with weak target-disease signal), acetazolamide → idiopathic intracranial hypertension (disease absent from OT entirely), verapamil → migraine (OT migraine associates only CACNA1A from the L-type channel family, not CACNA1C/D), anakinra → cryopyrin-associated periodic syndrome (IL1R1-CAPS subtype curation gap), and eculizumab → atypical haemolytic uremic syndrome (complement-gene subtype curation gap). The failures are not pipeline-logic errors but resource-coverage limits, which lift automatically with each OT release.

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


### Literature and patent prior separates crowded from underexplored frontiers

The literature pass ran across the top 5,000 ranked hypotheses, querying PubMed, Europe PMC, ClinicalTrials.gov, and the Lens.org patent backend in parallel. Of those 5,000 hypotheses, 1,883 (37.7%) carried PubMed evidence, 1,872 (37.4%) Europe PMC evidence, 1,177 (23.5%) at least one registered clinical trial, and 96 (1.9%) at least one Lens patent filing matched to the substance-disease pair; the four signals were log-saturated and combined into `investigation_prior ∈ [0, 1]` which damps opportunity by `(1 − 0.5 × investigation_prior)`. The patent prior carried real signal: the most-patented hypotheses in the top 5,000 were satralizumab → rheumatoid arthritis (125 filings, `investigation_prior` 0.99, demoted to rank 1,547), tralokinumab and lebrikizumab → psoriasis (20 filings each, prior 0.88), and the KRAS-G12C inhibitors adagrasib and sotorasib → cardiofaciocutaneous syndrome (18 filings each, prior 0.48) — separating mechanisms with material concurrent industry activity from genuinely under-investigated ones. The combined literature-and-patent prior carries through to the per-hit PDF briefs, where it is rendered as an explicit "investigation density" line so a reviewer can see at a glance whether the screen's hit sits in an active or quiet competitive space.

### Severity-aware damping eliminates the "receptor LoF + agonist" trap

Direct-target mechanism scoring produces a recurring failure mode in monogenic disorders: when the disease is caused by loss-of-function of a receptor, an agonist for that receptor cannot rescue the broken protein. Without correction, an early version of the screen had the insulin formulations clustered at ranks 11–30 of the full output paired with "hyperinsulinism due to INSR deficiency", and G-CSF analogues paired with "autosomal recessive severe congenital neutropenia due to CSF3R deficiency" — both biologically futile.

We implemented a narrow, conservative heuristic: when a gene named in the disease is also a direct target of the drug (`disease_gene_match`), the disease name contains severity language ("deficiency", "complete absence", "Donohue syndrome", "Rabson-Mendenhall"^26^), AND the drug's action class is agonistic, the hypothesis is flagged `severe_loF_agonist` and its opportunity is reduced by 70%. The cluster (32 hypotheses) was pushed below rank 100 — effectively out of the screen (**Fig. 3a**) — while adjacent partial-LoF rescuable cases were preserved: thiazolidinediones for PPARG-related familial partial lipodystrophy (FPLD3), where TZD agonism rescues residual receptor activity^27^, remained at rank 30; calcimimetics for familial hypocalciuric hypercalcemia 1, where partial CaSR function persists^28^, remained at ranks 19–20. The discrimination rests entirely on the disease name's severity language — INSR-deficiency hyperinsulinism implies a severe-LoF state, "partial lipodystrophy" or "hypocalciuric" does not.

### Filtering noisy direction signal from model-organism sources

A second class of failure mode arose from how Open Targets aggregates direction-of-effect evidence. Of approximately 1.02 million direction-informative rows in the OT evidence parquet, ~1.14 million originate from IMPC (mouse-knockout phenotyping via PhenoDigm), versus only ~166,000 from the human-genetics direction sources (ot_genetics_portal) and ~282,000 from ClinVar (eva) (**Fig. 3b**). IMPC's `variantEffect = LoF / directionOnTrait = risk` calls — appropriate for "what happens when you delete the mouse gene" — frequently flip the inferred therapeutic direction for human disorders whose mechanism differs from a complete knockout. For EPOR + primary familial polycythemia (a *gain*-of-function truncation that removes the receptor's negative-regulatory domain^29^), IMPC contributed 16 LoF/risk rows; the aggregated direction signal recommended an EPOR *agonist* — the opposite of the correct therapy. After filtering the direction adapter to seven curated human-genetics sources (`ot_genetics_portal`, `gene_burden`, `eva`, `gene2phenotype`, `clingen`, `genomics_england`, `orphanet`), the EPO → polycythemia cluster fell 275 ranks (22–30 → 297–301), and the direction CSV more broadly contracted from ~1.02M rows to 113,339 with a balanced direction split (64,294 +1 vs 49,045 −1, versus the previous ~95% skew toward +1).

![**Engineering choices that mattered.** **a**, Severity heuristic relocates all 32 'receptor LoF + agonist' hypotheses from their natural ranks (most of them clustering at ranks 10–30) to ranks below 100 (most below 30,000). The 0.70 damping factor was chosen to ensure no flagged hypothesis surfaces in a typical fund-screen review window (top 100). **b**, Direction-of-effect signal composition before and after the IMPC-source filter. IMPC mouse-knockout evidence (red) dominates the raw OT direction signal by an order of magnitude but is inappropriate for inferring human therapeutic direction in monogenic disorders; restricting to curated human-genetics sources (green) contracts the CSV from ~1.02M rows to 113k and rebalances the +1/−1 direction split from ~95% positive to ~57% positive.](figures/fig3_engineering_choices.png){ width=100% }


### Operational outputs

The screen produces a 36-column CSV that is hard for non-technical stakeholders to act on directly. We therefore added a one-page A4 PDF generator that, per hypothesis, weaves the mechanism evidence into prose ("Drug-target convergence on the disease scores 0.93; therapeutic direction is aligned with Open Targets genetic evidence; lead target is well-expressed in the disease-relevant tissue (brain); orthologous-gene model-organism evidence is strong"), translates the Orange Book columns into IP-runway commentary ("loss-of-exclusivity year 2037, with no generic competition — meaningful IP runway for a co-development conversation"), expresses market opportunity with explicit prevalence-source attribution ("approximately 16,650 US patients — below the 200,000 FDA Orphan Drug Act threshold; orphan exclusivity and premium pricing materially change the deal math"), and concludes with a proposed collaboration model under an investigator-initiated trial framework with the substance's IP holder. The brief is intended as a starting artefact for outreach to a clinical collaborator: it makes the mechanism case explicit, contextualises the regulatory and economic opportunity, and sketches a tractable next step.

To make the full ranked output explorable without writing Python, we ship a public Streamlit dashboard alongside the pipeline (Supplementary Fig. 2). The dashboard exposes the same 176,272-row output through three tabs: a *Browse* view (sidebar filters on drug / disease / target / opportunity / mech-support / orphan / IP runway / Lens patent activity; an opportunity-vs-mech scatter with click-through to a per-hit detail panel that renders the substance's 2D chemical structure from ChEMBL, a live PubMed link, and the bundled PDF brief); a *History* view (timeline of the 54 historic repurposings used to validate the screen, plus a disease-disease network where nodes are indications and edges are the drugs that bridge them — making the hub structure of historic repurposings, with hypertension, epilepsy, rheumatoid arthritis, and schizophrenia as the most fertile sources, immediately visible); and a *FAQ* view (one definition per computed value with the actual formula, designed so a non-technical collaborator can read a row without a tutorial). The dashboard is read-only, ships with a frozen parquet snapshot of the latest run, and is publicly available at **https://reprise.streamlit.app** (source mirror at https://github.com/wvcrieki42/reprise).

## Discussion

The technical content of each enrichment layer in isolation is not novel: novelty subtraction via ontology hops, directionality inference from genetic evidence, tissue-expression gating, pathway-overlap analysis, and PubMed-based investigation priors are all individually established. The contribution of this work is the *integration* — combining these layers into a single pipeline with measured precision (89% recovery of known repurposing wins, n = 54), with the regulatory and economic signals (IP runway, prevalence, generic status, patent activity) that determine whether a candidate is actionable, and with the operational artefact (PDF outreach briefs) that determines whether the work leaves the screen.

The validation framework is also worth emphasising. Many repurposing screens are presented with anecdotal example hits; here we measure precision against an explicit ground-truth set of 54 FDA / EMA approvals and well-established off-label uses, characterise every failure mode (all six misses are OT data-coverage gaps, each traceable to a specific missing target-disease association — e.g. verapamil → migraine where OT associates only CACNA1A from the L-type channel family, not the CACNA1C/D channels verapamil actually hits clinically), and ship the measurement script alongside the pipeline so any user can re-verify or extend. The metformin/PCOS surprising-hit case is, in our view, a small piece of evidence that the screen is doing genuine mechanistic inference rather than retrieval — it surfaces the proximal pharmacological mechanism rather than the post-hoc literature attribution.

The combination-therapy result deserves note. Re-deriving the BRAF + MEK doublet for cardiofaciocutaneous syndrome from first principles — without any oncology-developmental cross-reference programmed in — is a particularly clean form of internal validation. The screen's noisy-OR formulation across pathway-bridging companion targets recovers an established clinical strategy in an adjacent indication.

Limitations are real. Open Targets' disease-association coverage is the precision ceiling, and our six backtest misses correspond exactly to under-curated rare or peripheral indications (each traceable to either a missing canonical mechanism target in OT or to a disease subtype not yet propagated into the target-disease join). The severity heuristic relies on disease-name keywords rather than structured mechanism-of-disease annotations from, e.g., ClinGen; expanding to structured severity will require curated input. The literature pass currently covers the top 5,000 hypotheses; extending it across the full 176,272-row output would require either a multi-day cache build or a more parallel scheduler than the current rate-limiter permits. The Orange Book is small-molecule only; the FDA Purple Book equivalent for biologics is an obvious extension.

The pipeline is engineered for incremental addition. Natural next layers include cell-line / DepMap evidence for cancer hypotheses^33^, structured drug-drug interaction filtering for combination-therapy safety^34^, ClinGen-derived structured severity annotations to replace the current disease-name keyword heuristic, an FDA Purple Book equivalent for biologics IP runway (which currently fall through the Orange Book layer), and a nightly refresh pipeline that compares each new upstream release for regressions before promoting the run output into the public dashboard.

The full pipeline source, configuration, curated CSVs, validation YAML, and this manuscript are publicly available at https://github.ugent.be/wvcrieki/repurposed.

## Methods (condensed)

**Data sources.** ChEMBL 37^14^ (approved drugs, mechanism of action, indications, molecule_hierarchy), Open Targets 24.06^12^ (target-disease associations, evidence, diseases, baseline expression, targets parquet with integrated Reactome^30^ pathways), STRING v12.0^31^ (protein-protein interaction edges), EFO 2024-06^32^ (disease ontology), FDA Orange Book (May 2026 release), Orphanet (en_product9_prev.xml, CC-BY-4.0)^18^, PubMed E-utilities, Europe PMC REST, ClinicalTrials.gov v2 API, Lens.org Patent API.

**Mechanistic scoring.** Per (drug, disease), `mech_support = 1 − ∏ᵢ (1 − contribᵢ)` (noisy-OR^19^) where each contribution is `target_weight × assoc_score` over the drug's direct targets (weight 1.0) optionally augmented by STRING neighbours (weight 0.5 × STRING_score). Open Targets `assoc_score` filtered at >= 0.1.

**Novelty.** Ontology distance from the drug's nearest known indication via EFO graph^32^ walked outward to radius 1. ChEMBL drug indications normalised colon→underscore (`MONDO:0005148 → MONDO_0005148`) to align with OT format. Target-class indication rollup: drug A inherits drug B's indications iff A's frozenset of (direct target, action_type) pairs is a subset of B's; STRING neighbours excluded from class keys.

**Direction.** Filtered Open Targets evidence to `ot_genetics_portal`, `gene_burden`, `eva`, `gene2phenotype`, `clingen`, `genomics_england`, `orphanet` only (IMPC and other non-human-genetics sources excluded). Per (target, disease), `variantEffect ∈ {LoF, GoF}` combined with `directionOnTrait ∈ {risk, protective}` gives therapeutic direction in {+1, −1} per the OT documented inference^12^.

**Tissue.** Per (target, disease), score is max over the disease's relevant tissues of (`disease_relevance × target_expression`). Bucketed into expressed (>= 0.25), low (linear ramp), absent (measured zero, factor 0.3), unknown (no measurement, factor 0.7). Multiplicative.

**Phylogenetic boost.** From OT evidence `datatypeId == animal_model` (default IMPC^16^), per (drug, disease) the max `phylo_score` across the drug's targets. `phylo_factor = 1 + 0.5 × phylo_score`. Asymmetric: presence boosts, absence does not penalise.

**Pathway boost.** Reactome^30^ pathway co-membership per (drug, disease), restricted to pathways with ≤ 80 genes (specificity filter). Only drug_target ≠ disease_target bridges credited (`indirect_only = true`, to avoid double-credit of direct hits). `pathway_factor = 1 + 0.3 × log1p(n_overlap) / log1p(5)`.

**Opportunity score.** `mech^1 × novelty^1 × direction × tissue × phylo × pathway / sqrt(n_drug_targets)`.

**Substance grouping.** ChEMBL `molecule_hierarchy.parent_molregno` rolls formulations up to active ingredient.

**Severity damping.** When `disease_gene_match` populated AND disease name contains severity language AND drug action is agonistic, opportunity × 0.3, re-rank.

**Literature pass.** Per top-N pair (default 5,000), query PubMed/Europe PMC/ClinicalTrials.gov/Lens with synonym-expanded queries. Counts log-saturated and combined to `investigation_prior ∈ [0, 1]`. Opportunity damped by `(1 − 0.5 × investigation_prior)`.

**Market sizing.** Curated CSV (122 major US diseases from CDC/SEER/foundation estimates) overlaid by Orphanet structured prevalence (~5,100 rare diseases bridged to OT MONDO via dbXRefs).

**IP signals.** FDA Orange Book parsed per active ingredient: max patent expiry, max exclusivity expiry, generic availability via any ANDA filing.

**Backtest validation.** 54 curated cases (FDA / EMA approvals and well-established off-label uses across oncology, cardiometabolic, neuro / psychiatry, immunology, and rare disease) with expected drug, repurposed indication, and mechanism targets. Per case, all ChEMBL IDs matching drug name (parent + salts) pooled; mech_support computed; HIT at >= 0.3. The 95% Clopper–Pearson confidence interval on the hit rate is 77–96%.

**Output.** 36-column CSV plus one-page A4 PDFs (ReportLab) per top hypothesis that translate the row into a mechanism-rationale brief with IP and market context.

**Dashboard.** Public Streamlit application (Browse / History / FAQ tabs) bundled with a frozen parquet snapshot of the full-screen output and the 30 highest-ranked PDF briefs. The detail panel fetches 2D structures live from ChEMBL's image endpoint, cached 24 h. The history-tab network is computed via `networkx` spring layout over (original_indication, repurposed_indication) edges drawn from the curated validation YAML.

**Compute.** Full run with all enrichment layers and literature pass: approximately 12 minutes on M-series Mac (16 GB RAM, ~50 GB free disk). DuckDB engine for the streaming target-disease join^35^.

## Tables

**Table 1.** Representative subset (19 of 54) of the curated backtest set of known successful drug-repurposing wins, with the original indication, the repurposed indication scored by the screen, the expected mechanism targets, the `mech_support` value computed by the screen, and the HIT / MISS / DATA-GAP status at threshold 0.30. The remaining 35 cases (oncology, cardiometabolic, neuro/psych, immunology, and rare-disease extensions, full list in `data/curated/repurposing_validation.yaml`) follow the same pattern: 32 HIT, 3 OT gap.

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
| | | | **Subset HIT rate** | **16/19** | **84%** |
| | | | **Full-set HIT rate** | **48/54** | **89%** |

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

## Supplementary materials

**Supplementary Fig. 1 — Example one-page collaboration brief** (`supplementary/supplementary_fig1_example_brief.pdf`). Auto-generated brief for the screen's top-ranked hypothesis, bromazepam → developmental and epileptic encephalopathy, exported by the `build_deal_memos.py --no-kol` profile used for blinded outreach. The brief weaves the mechanism evidence, the IP-runway signal, the orphan-prevalence opportunity, and a proposed investigator-initiated trial framework into a single A4 page. The same template generates briefs for any row of the 36-column output CSV.

**Supplementary Fig. 2 — Public Streamlit dashboard** (live at **https://reprise.streamlit.app**, source at https://github.com/wvcrieki42/reprise and `dashboard/` in the primary repository). Interactive companion to the manuscript that exposes the full 176,272-row ranked output through three tabs (Browse / History / FAQ), with sidebar filters, an opportunity-vs-mech scatter, per-hit detail panels rendering 2D structures from ChEMBL, live PubMed search links, and the bundled PDF briefs. The History tab carries the timeline of the 54-case validation set and the disease-disease network. The dashboard is hosted on Streamlit Community Cloud and bundles a frozen parquet snapshot of the screen's full output.

## Acknowledgements

This work was performed at BioBix (Department of Mathematical Modelling, Statistics and Bioinformatics, Ghent University, Belgium). We thank the Open Targets, ChEMBL, Reactome, EFO, STRING, IMPC, Orphanet, FDA Orange Book, PubMed, Europe PMC, ClinicalTrials.gov, and Lens.org teams for maintaining and openly releasing the data resources on which this work depends.

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
16. Muñoz-Fuentes, V. *et al.* The International Mouse Phenotyping Consortium (IMPC): a functional catalogue of the mammalian genome that informs conservation. *Conserv. Genet.* **19**, 995–1005 (2018).
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
30. Gillespie, M. *et al.* The Reactome pathway knowledgebase 2022. *Nucleic Acids Res.* **50**, D687–D692 (2022).
31. Szklarczyk, D. *et al.* The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Res.* **51**, D638–D646 (2023).
32. Malone, J. *et al.* Modeling sample variables with an Experimental Factor Ontology. *Bioinformatics* **26**, 1112–1118 (2010).
33. Tsherniak, A. *et al.* Defining a cancer dependency map. *Cell* **170**, 564–576 (2017).
34. Wishart, D. S. *et al.* DrugBank 5.0: a major update to the DrugBank database for 2018. *Nucleic Acids Res.* **46**, D1074–D1082 (2018).
35. Raasveldt, M. & Mühleisen, H. DuckDB: an embeddable analytical database. In *Proc. ACM SIGMOD International Conference on Management of Data* 1981–1984 (2019).
