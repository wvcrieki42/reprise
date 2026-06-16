---
geometry: margin=2.4cm
fontsize: 11pt
mainfont: "Times New Roman"
---

\noindent\textbf{Prof. Dr. Wim Van Criekinge}\
BioBix — Department of Mathematical Modelling, Statistics and Bioinformatics\
Faculty of Bioscience Engineering, Ghent University\
Coupure Links 653, B-9000 Ghent, Belgium\
Wim.VanCriekinge@UGent.be

\vspace{0.6cm}

\noindent The Editor\
\textit{Nature Biotechnology}

\vspace{0.4cm}

\noindent Dear Editor,

I would like to submit the enclosed Article, **"REPRISE: integrated mechanism, regulatory and prevalence scoring for drug repurposing"**, for consideration in *Nature Biotechnology*. The manuscript describes a public computational platform that scores all approved drugs against the disease ontology in a single pass, with measured precision and explicit translational context, and it ships with a live, openly accessible dashboard at **https://reprise.streamlit.app**.

The contribution we believe is novel and well-suited to your readership:

- **Integrated screening at scale.** REPRISE evaluates **176,272 drug–disease hypotheses** (3,996 approved drugs × 28,198 ontology diseases) in approximately **12 minutes on a laptop**, layering noisy-OR mechanism scoring with directionality from genetic evidence, tissue expression, Reactome pathway co-membership, FDA Orange Book IP runway, Orphanet rare-disease prevalence, and a four-source literature-and-patent prior (PubMed, Europe PMC, ClinicalTrials.gov, Lens.org) computed across **every ranked hypothesis** — not just the top *N*. Of the 176,272 hypotheses, 83,042 carry PubMed evidence, 90,123 Europe PMC evidence, 61,721 a registered clinical trial, and 3,008 a Lens patent filing.

- **Measured precision, not anecdote.** Against a curated set of **54 historic FDA / EMA repurposing successes**, REPRISE recovers **48 (89 %, 95 % Clopper–Pearson CI 77–96 %)** at a mechanism-support threshold of 0.30. Every one of the six misses traces to a specific Open Targets coverage gap, not a pipeline-logic error — a transparent, reproducible failure analysis that we believe sets a new bar for repurposing-screen validation.

- **A clinically meaningful** *de novo* **finding.** Without any prior cross-reference between oncology and developmental medicine, the screen rediscovers the **BRAF + MEK doublet** — currently standard-of-care for *BRAF*-V600E metastatic melanoma — as a candidate combination therapy for **cardiofaciocutaneous syndrome**, a rare developmental RASopathy. The screen also surfaces metformin's mitochondrial complex-I mechanism in polycystic ovary syndrome, recovering the proximal pharmacological pathway rather than the AMPK route the literature most often invokes.

- **Built for action.** Each top hypothesis ships as a **one-page PDF brief** that translates the row of numbers into a narrative for clinical-collaborator outreach, and the live Streamlit dashboard exposes the full ranked output with chemical structures, live PubMed links, IP signals, and a disease-disease network making the historical hub structure of repurposing (hypertension, epilepsy, type 2 diabetes, rheumatoid arthritis and schizophrenia together supplied >40 % of curated wins) immediately visible.

Beyond academic interest, REPRISE was designed to be operationally useful: BioBix at Ghent University and partners in venture capital and clinical research currently use it as an early-stage screening layer for repurposing programmes. We believe the integration of measured mechanistic precision with explicit regulatory and economic context — and the open public dashboard — fits the *Nature Biotechnology* mandate of rigorous, broadly applicable biotechnology platforms.

This manuscript has not been published and is not under consideration at any other journal. The source code, configuration, curated validation set, run output and this manuscript are publicly available at **https://github.ugent.be/wvcrieki/repurposed**. I have no competing interests to declare.

I would be happy to suggest expert reviewers if helpful.

Thank you for considering the manuscript.

\vspace{0.4cm}

\noindent Yours sincerely,

\vspace{0.8cm}

\noindent **Prof. Dr. Wim Van Criekinge**\
Full Professor, Department of Mathematical Modelling, Statistics and Bioinformatics\
Ghent University
