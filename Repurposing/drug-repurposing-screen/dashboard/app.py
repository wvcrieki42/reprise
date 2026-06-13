"""REPRISE -- public companion dashboard.

A Streamlit app that browses the 176,272 ranked drug-repurposing hypotheses
produced by the REPRISE screen. Designed as a paper companion: read-only,
no secrets, deployable to Streamlit Community Cloud or HF Spaces in one click.

Run locally:
    streamlit run dashboard/app.py

Data file: dashboard/data/repurposing_hypotheses.parquet (snapshot bundled
with the repo). To refresh, copy the latest output_full/...parquet over it.
"""
from __future__ import annotations

import base64
import re
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import yaml

DATA = Path(__file__).parent / "data" / "repurposing_hypotheses.parquet"
VALIDATION_YAML = Path(__file__).parent / "data" / "repurposing_validation.yaml"
BRIEFS_DIR = Path(__file__).parent / "static" / "briefs"
# Streamlit serves files from dashboard/static/ at /app/static/<file> when
# enableStaticServing = true (see .streamlit/config.toml). The relative URL
# below works both locally and on Streamlit Community Cloud.
BRIEFS_URL_PREFIX = "/app/static/briefs/"
PAPER_URL = "https://github.ugent.be/wvcrieki/repurposed"


def _safe_filename(s: str) -> str:
    """Mirror scripts/build_deal_memos.py::_safe_filename so the lookup matches."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s.strip())[:80] or "memo"


def brief_path_for(substance: str, disease: str) -> Path | None:
    """Return the path to a bundled brief PDF for this (substance, disease)
    if one exists in dashboard/briefs/, else None."""
    if not substance or not disease:
        return None
    p = BRIEFS_DIR / f"{_safe_filename(str(substance))}__{_safe_filename(str(disease))}__us.pdf"
    return p if p.exists() else None


@st.cache_data(show_spinner=False)
def brief_data_uri(brief_path_str: str) -> str:
    """Encode a brief PDF as a data: URI so LinkColumn opens it without
    depending on Streamlit's static-file serving."""
    p = Path(brief_path_str)
    if not p.exists():
        return ""
    b64 = base64.b64encode(p.read_bytes()).decode("ascii")
    return f"data:application/pdf;base64,{b64}"


CHEMBL_IMG_URL = "https://www.ebi.ac.uk/chembl/api/data/image/{chembl_id}.svg"


@st.cache_data(show_spinner=False)
def load_repurposing_history() -> pd.DataFrame:
    """Load the 54-case curated YAML as a DataFrame for the History tab."""
    if not VALIDATION_YAML.exists():
        return pd.DataFrame()
    cases = yaml.safe_load(VALIDATION_YAML.read_text())["cases"]
    rows = []
    for c in cases:
        rows.append({
            "drug": c["drug"].title(),
            "original": str(c.get("original_indication", "")).strip(),
            "repurposed": str(c.get("repurposed_disease_name", "")).strip(),
            "year": c.get("year_repurposed"),
            "targets": ", ".join(c.get("mechanism_targets", []) or []),
            "notes": c.get("notes", ""),
        })
    df = pd.DataFrame(rows)
    df = df[df["year"].notna()].copy()
    df["year"] = df["year"].astype(int)
    return df


@st.cache_data(show_spinner=False, ttl=24 * 3600)
def fetch_chembl_structure(chembl_id: str) -> str | None:
    """Fetch a 2D structure SVG from ChEMBL's public image endpoint.
    Returns the raw SVG text on success, None for biologics / unknown IDs
    (ChEMBL returns 400 for entries with no 2D structure)."""
    if not chembl_id or not str(chembl_id).startswith("CHEMBL"):
        return None
    import urllib.request
    import urllib.error
    try:
        req = urllib.request.Request(
            CHEMBL_IMG_URL.format(chembl_id=chembl_id),
            headers={"User-Agent": "REPRISE-dashboard/1.0"})
        with urllib.request.urlopen(req, timeout=8) as resp:
            if resp.status != 200:
                return None
            ctype = resp.headers.get("Content-Type", "")
            if "svg" not in ctype:
                return None
            return resp.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ValueError):
        return None

st.set_page_config(
    page_title="REPRISE -- drug-repurposing screen",
    page_icon=":pill:",
    layout="wide",
)


# ----------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------
PUBMED_SEARCH_URL = "https://pubmed.ncbi.nlm.nih.gov/?term="


def _pubmed_query_url(target: str, disease: str, count: int) -> str:
    """Build a PubMed search URL for this (target, disease) pair. The count is
    embedded as a fragment so the LinkColumn display_text regex can render it
    as the clickable label without changing the URL semantics."""
    import urllib.parse
    if not target or not disease or not count:
        return ""
    term = f'"{target}"[tiab] AND "{disease}"[tiab]'
    return f"{PUBMED_SEARCH_URL}{urllib.parse.quote_plus(term)}#n={int(count)}"


@st.cache_data(show_spinner=False)
def load_data() -> pd.DataFrame:
    df = pd.read_parquet(DATA)
    for col in ("opportunity", "mechanistic_support", "novelty",
                "investigation_prior", "pathway_score"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    # Pre-compute a row-level link to the bundled PDF brief when one exists.
    # Embed as a base64 data: URI -- the browser opens it directly, no
    # dependency on Streamlit's static-file serving (which is finicky to
    # configure across local + Cloud + reverse-proxy deploys).
    name_col = "substance_name" if "substance_name" in df.columns else "drug_name"
    brief_urls = []
    for s, d in zip(df[name_col].fillna(""), df["disease_name"].fillna("")):
        p = brief_path_for(s, d)
        brief_urls.append(brief_data_uri(str(p)) if p is not None else "")
    df["brief_url"] = brief_urls
    # PubMed search URL per row -- non-empty only when count > 0 so the
    # LinkColumn renders blank for 0/NA rows and a clickable count for hits.
    if {"pubmed_count", "lead_target", "disease_name"} <= set(df.columns):
        df["pubmed_link"] = [
            _pubmed_query_url(t, d, c) if pd.notna(c) and int(c) > 0 else ""
            for t, d, c in zip(df["lead_target"].fillna(""),
                               df["disease_name"].fillna(""),
                               df["pubmed_count"].fillna(0))
        ]
    return df


df = load_data()


# ----------------------------------------------------------------------
# Header
# ----------------------------------------------------------------------
st.title("REPRISE")
st.caption(
    "**Repurposing Engine for Pathway-Resolved Indication Scoring and Evidence** -- "
    "a mechanism-driven screen across all approved drugs and disease ontology."
)

n_total = len(df)
n_orphan = int(df["is_orphan"].fillna(False).astype(bool).sum()) if "is_orphan" in df.columns else 0
n_drugs = df["substance_chembl_id"].nunique() if "substance_chembl_id" in df.columns else 0
n_disease = df["efo_id"].nunique() if "efo_id" in df.columns else 0

c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Hypotheses", f"{n_total:,}")
c2.metric("Substances", f"{n_drugs:,}")
c3.metric("Diseases", f"{n_disease:,}")
c4.metric("Orphan-flagged", f"{n_orphan:,}")
c5.metric("Backtest precision", "89%",
          help="48/54 known repurposings recovered, 95% CI 77-96%")


tab_browse, tab_history, tab_faq = st.tabs([
    "Browse hypotheses",
    "History & how REPRISE adds value",
    "FAQ: how to read this dashboard",
])


# ======================================================================
# Tab 1 -- Browse
# ======================================================================
with tab_browse:
    # ------------------------------------------------------------------
    # Sidebar filters (shared, but only meaningful while Browse is active)
    # ------------------------------------------------------------------
    st.sidebar.header("Filters")

    drug_q = st.sidebar.text_input(
        "Active ingredient contains",
        help="Substance name after rolling up ChEMBL formulations to the active ingredient."
    ).strip().upper()
    disease_q = st.sidebar.text_input("Disease name contains").strip().lower()
    target_q = st.sidebar.text_input("Lead target contains").strip().upper()

    opp_min, opp_max = float(df["opportunity"].min()), float(df["opportunity"].max())
    opp_lo, opp_hi = st.sidebar.slider(
        "Opportunity score",
        min_value=round(opp_min, 2), max_value=round(opp_max, 2),
        value=(round(opp_min, 2), round(opp_max, 2)), step=0.05,
    )

    if "mechanistic_support" in df.columns:
        mech_lo, mech_hi = st.sidebar.slider(
            "Mechanism support",
            min_value=0.0, max_value=1.0,
            value=(0.0, 1.0), step=0.05,
        )
    else:
        mech_lo, mech_hi = 0.0, 1.0

    orphan_only = st.sidebar.checkbox("Orphan indications only", value=False)
    has_ip = st.sidebar.checkbox("Has Orange Book IP runway", value=False)
    has_patents = st.sidebar.checkbox("Has Lens patent activity", value=False)
    exclude_severity = st.sidebar.checkbox(
        "Hide severity-flagged (LoF + agonist)", value=True,
        help="Suppresses the curated 'receptor LoF + agonist' false-positive cluster.")

    st.sidebar.divider()
    st.sidebar.caption(
        f"Source code & paper: [github.ugent.be/wvcrieki/repurposed]({PAPER_URL})"
    )

    # ------------------------------------------------------------------
    # Filter application
    # ------------------------------------------------------------------
    mask = pd.Series(True, index=df.index)
    if drug_q:
        name_col = "substance_name" if "substance_name" in df.columns else "drug_name"
        mask &= df[name_col].fillna("").str.upper().str.contains(drug_q, regex=False)
    if disease_q:
        mask &= df["disease_name"].fillna("").str.lower().str.contains(disease_q, regex=False)
    if target_q and "lead_target" in df.columns:
        mask &= df["lead_target"].fillna("").str.upper().str.contains(target_q, regex=False)
    mask &= df["opportunity"].between(opp_lo, opp_hi)
    if "mechanistic_support" in df.columns:
        mask &= df["mechanistic_support"].fillna(0).between(mech_lo, mech_hi)
    if orphan_only and "is_orphan" in df.columns:
        mask &= df["is_orphan"].fillna(False).astype(bool)
    if has_ip and "latest_patent_year" in df.columns:
        mask &= df["latest_patent_year"].notna()
    if has_patents and "patent_count" in df.columns:
        mask &= df["patent_count"].fillna(0).gt(0)
    if exclude_severity and "severity_concern" in df.columns:
        mask &= ~df["severity_concern"].fillna("").astype(str).str.contains(
            "severe_loF_agonist", case=False, regex=False)

    view = df[mask].sort_values("rank")

    # ------------------------------------------------------------------
    # Results table
    # ------------------------------------------------------------------
    st.subheader(f"Ranked hypotheses ({len(view):,} match{'es' if len(view) != 1 else ''})")

    if view.empty:
        st.info("No hypotheses match these filters. Loosen them in the sidebar.")
        st.stop()

    # ------------------------------------------------------------------
    # Scatter -- opportunity vs mech, click any point to drill in
    # ------------------------------------------------------------------
    scatter_n = min(2000, len(view))
    scatter_view = view.head(scatter_n).reset_index().rename(columns={"index": "_df_idx"})
    fig = px.scatter(
        scatter_view,
        x="opportunity",
        y="mechanistic_support",
        color=scatter_view["is_orphan"].fillna(False).astype(bool).map(
            {True: "orphan", False: "non-orphan"}),
        color_discrete_map={"orphan": "#D55E00", "non-orphan": "#0072B2"},
        hover_data={
            "substance_name": True,
            "disease_name": True,
            "lead_target": True,
            "rank": True,
            "opportunity": ":.3f",
            "mechanistic_support": ":.3f",
            "_df_idx": False,
        },
        labels={"opportunity": "Opportunity (composite)",
                "mechanistic_support": "Mechanism support (noisy-OR)",
                "color": ""},
        title=f"Opportunity vs mechanism support  (top {scatter_n} of "
              f"{len(view):,} matches; click any point)",
    )
    fig.update_traces(marker=dict(size=7, opacity=0.75, line=dict(width=0)))
    fig.update_layout(height=460, margin=dict(l=10, r=10, t=50, b=10),
                      legend=dict(orientation="h", yanchor="bottom", y=1.02,
                                  xanchor="right", x=1))
    scatter_event = st.plotly_chart(
        fig, use_container_width=True, key="scatter",
        on_select="rerun", selection_mode="points",
    )

    display_cols = [
        c for c in [
            "rank", "brief_url", "substance_name", "disease_name", "lead_target",
            "opportunity", "mechanistic_support", "novelty", "direction_status",
            "tissue_status", "is_orphan", "us_patients",
            "latest_patent_year", "has_generic",
            "pubmed_link", "trial_count", "patent_count", "investigation_prior",
            "combo_partner_1_name", "combo_partner_1_synergy",
        ] if c in view.columns
    ]

    # --- Explicit fallback selector for the per-hit detail panel ---
    # If the user's environment doesn't fire dataframe row-click events
    # (some Streamlit versions / browser configs) this dropdown still works.
    rows_with_brief = view[view["brief_url"] != ""].head(30) if "brief_url" in view.columns else view.head(0)
    if not rows_with_brief.empty:
        st.markdown("### Inspect a top hypothesis")
        st.caption(
            "Pick from the dropdown OR click any row / scatter point below. "
            "The detail panel renders the chemical structure and embeds the "
            "PDF brief inline."
        )
        labels = [
            f"#{int(r['rank']):>3}  {r['substance_name']} -> {r['disease_name']}"
            for _, r in rows_with_brief.iterrows()
        ]
        index_map = list(rows_with_brief.index)
        choice = st.selectbox(
            "Hypothesis", options=["(none -- click row or scatter point below)"] + labels,
            index=0, key="hyp_picker", label_visibility="collapsed",
        )
        if choice != "(none -- click row or scatter point below)":
            picked_index = index_map[labels.index(choice)]
            st.session_state["picker_index"] = int(picked_index)
        else:
            st.session_state["picker_index"] = None
    else:
        st.session_state["picker_index"] = None

    event = st.dataframe(
        view[display_cols].head(2000),
        use_container_width=True,
        hide_index=True,
        column_config={
            "rank": st.column_config.NumberColumn("rank", format="%d",
                help="Position in the screen's global ranking by opportunity."),
            "substance_name": st.column_config.TextColumn("Active ingredient",
                help="ChEMBL active-ingredient name after rolling salts and "
                     "formulations up to the molecule_hierarchy parent."),
            "disease_name": st.column_config.TextColumn("Disease",
                help="Open Targets / EFO canonical disease name."),
            "lead_target": st.column_config.TextColumn("Lead target",
                help="Highest-contributing drug target for this hypothesis."),
            "opportunity": st.column_config.NumberColumn("Opportunity", format="%.3f",
                help="Composite score: mech x novelty x direction x tissue x phylo x "
                     "pathway / sqrt(n_targets). Higher = more attractive."),
            "mechanistic_support": st.column_config.NumberColumn("Mech", format="%.3f",
                help="Noisy-OR aggregation of drug-target x target-disease evidence. "
                     "0 = no mechanism, 1 = saturated coverage."),
            "novelty": st.column_config.NumberColumn("Novelty", format="%.2f",
                help="1 if the indication is not in the drug's ChEMBL on-label list "
                     "and not in any same-target same-action class member's list."),
            "direction_status": st.column_config.TextColumn("Direction",
                help="aligned / opposing / unknown -- does the drug push the target "
                     "the way OT genetics says is therapeutic?"),
            "tissue_status": st.column_config.TextColumn("Tissue",
                help="expressed / low / absent / unknown -- is the target present in "
                     "the disease-relevant tissue?"),
            "is_orphan": st.column_config.CheckboxColumn("Orphan",
                help="< 200,000 US patients (Orphan Drug Act threshold)."),
            "us_patients": st.column_config.NumberColumn("US patients", format="%d"),
            "latest_patent_year": st.column_config.NumberColumn("Patent yr", format="%d",
                help="Latest FDA Orange Book patent expiry year for the active ingredient."),
            "has_generic": st.column_config.CheckboxColumn("Generic",
                help="True if any ANDA generic has been filed for this active ingredient."),
            "pubmed_link": st.column_config.LinkColumn(
                "PubMed",
                help="Number of PubMed papers co-mentioning the lead target and "
                     "the disease in title/abstract. Click to open the search "
                     "in a new tab.",
                display_text=r"#n=(\d+)$",
            ),
            "trial_count": st.column_config.NumberColumn("Trials", format="%d"),
            "patent_count": st.column_config.NumberColumn("Patents", format="%d",
                help="Lens.org patent filings matching the target + disease search."),
            "investigation_prior": st.column_config.NumberColumn("Inv. prior", format="%.2f",
                help="[0,1] -- log-saturated aggregation across PubMed / Europe PMC / "
                     "NCT / Lens. Higher = more crowded frontier; damps opportunity."),
            "combo_partner_1_name": st.column_config.TextColumn("Combo partner",
                help="Top combination companion drug -- bridges a disease target the "
                     "primary substance misses."),
            "combo_partner_1_synergy": st.column_config.NumberColumn("Synergy", format="%.3f",
                help="combo_mech_support - primary_mech_support under the noisy-OR."),
            "brief_url": st.column_config.LinkColumn(
                ":page_facing_up: Brief",
                help="One-page PDF brief: mechanism rationale, IP runway, "
                     "clinical opportunity, and proposed collaboration model. "
                     "Click to open in a new tab. Bundled for the top-30 "
                     "ranked hypotheses; empty for the rest.",
                display_text="open PDF",
                pinned=True,
                width="small",
            ),
        },
        on_select="rerun",
        selection_mode="single-row",
        key="hits_table",
    )

    st.caption("Showing the top 2,000 rows after filtering. Refine in the sidebar to "
               "drill in. Hover any column header for its definition. Click any row "
               "(or any scatter point above) for the full per-hypothesis brief.")

    # ------------------------------------------------------------------
    # Per-hit detail panel
    # ------------------------------------------------------------------
    def _fmt_int(v) -> str:
        try:
            return f"{int(float(v)):,}"
        except Exception:
            return "n/a"

    def render_brief(row: pd.Series) -> None:
        substance = str(row.get("substance_name") or row.get("drug_name") or "?")
        disease = str(row.get("disease_name") or "?")
        chembl_id = str(row.get("substance_chembl_id") or row.get("drug_id") or "")

        st.subheader(f"{substance}  ->  {disease}")
        if chembl_id.startswith("CHEMBL"):
            st.caption(
                f"ChEMBL ID: [{chembl_id}](https://www.ebi.ac.uk/chembl/"
                f"explore/compound/{chembl_id})"
            )

        # ----- Chemical structure (sandboxed iframe -- bypasses every Streamlit
        # HTML sanitiser, works across all versions and deploy modes) ----------
        if chembl_id.startswith("CHEMBL"):
            svg = fetch_chembl_structure(chembl_id)
            if svg:
                # Strip ChEMBL's fixed 500x500 dimensions and re-set so the
                # SVG fits the iframe and centers cleanly.
                svg_fit = re.sub(r"\swidth='[^']*'", " width='360'", svg, count=1)
                svg_fit = re.sub(r"\sheight='[^']*'", " height='360'", svg_fit, count=1)
                struct_html = (
                    "<!doctype html><html><head><meta charset='utf-8'>"
                    "<style>body{margin:0;display:flex;justify-content:center;"
                    "align-items:center;background:#fff;font-family:system-ui}"
                    ".box{border:1px solid #e0e0e0;border-radius:8px;padding:10px;"
                    "background:#fff;text-align:center}"
                    ".cap{font-size:11px;color:#666;margin-top:4px}</style></head>"
                    f"<body><div class='box'>{svg_fit}"
                    "<div class='cap'>2D structure (ChEMBL)</div></div></body></html>"
                )
                st.components.v1.html(struct_html, height=410, scrolling=False)
            else:
                st.caption(
                    ":dna: _2D structure not available -- biologic or other "
                    "non-small-molecule substance (mAb, peptide, gene therapy)._"
                )

        st.markdown("**Why this hypothesis is real**")
        bits = []
        mech = row.get("mechanistic_support")
        if pd.notna(mech):
            bits.append(f"Drug-target convergence on the disease scores **{float(mech):.2f}** "
                        "(noisy-OR over Open Targets target-disease associations).")
        nov = row.get("novelty_status") or row.get("novelty")
        if isinstance(nov, str):
            bits.append(f"Novelty status: **{nov}**.")
        direction = row.get("direction_status")
        if direction and direction != "unknown":
            bits.append(f"Therapeutic direction: **{direction}**.")
        tissue = row.get("tissue_status")
        if tissue and tissue != "unknown":
            bits.append(f"Tissue: lead target is **{tissue}** in the disease-relevant tissue.")
        phylo = row.get("phylo_score")
        if pd.notna(phylo) and float(phylo) > 0:
            bits.append(f"Phylogenetic evidence: **{float(phylo):.2f}** "
                        f"({_fmt_int(row.get('phylo_n_models'))} model-organism rows).")
        n_overlap = row.get("n_pathway_overlap")
        if pd.notna(n_overlap) and int(n_overlap) > 0:
            bits.append(f"Reactome pathway co-membership: **{int(n_overlap)}** shared specific pathway(s).")
        st.markdown("  \n".join(f"- {b}" for b in bits) if bits else "_no enrichment data_")

        st.markdown("**Why now (IP runway + literature & patent prior)**")
        ip_bits = []
        if pd.notna(row.get("latest_patent_year")):
            ip_bits.append(f"Latest Orange Book patent expiry: **{int(row['latest_patent_year'])}**.")
        if pd.notna(row.get("loe_year")):
            ip_bits.append(f"Loss-of-exclusivity year: **{int(row['loe_year'])}**.")
        if row.get("has_generic"):
            ip_bits.append("Generic competitor available (ANDA filed).")
        elif pd.notna(row.get("latest_patent_year")):
            ip_bits.append("No ANDA generic on file.")
        else:
            ip_bits.append("Outside the Orange Book (biologic or other) -- check Purple Book.")
        pm = row.get("pubmed_count")
        if pd.notna(pm):
            ip_bits.append(f"PubMed: {_fmt_int(pm)} co-mentions; ClinicalTrials.gov: "
                           f"{_fmt_int(row.get('trial_count'))} registered trials; "
                           f"Lens patents: {_fmt_int(row.get('patent_count'))}.")
        if pd.notna(row.get("investigation_prior")):
            ip_bits.append(f"Aggregated investigation prior: **{float(row['investigation_prior']):.2f}** "
                           "(higher = more crowded frontier).")
        st.markdown("  \n".join(f"- {b}" for b in ip_bits))

        st.markdown("**Clinical opportunity**")
        if pd.notna(row.get("us_patients")):
            pts = int(row["us_patients"])
            if row.get("is_orphan"):
                st.markdown(
                    f"- US population: approximately **{pts:,}** patients -- below the "
                    f"200,000 FDA Orphan Drug Act threshold. Orphan exclusivity + premium "
                    f"pricing changes the deal math materially."
                )
            else:
                st.markdown(f"- US population: approximately **{pts:,}** patients.")
            if row.get("market_source"):
                st.caption(f"source: {row['market_source']}")
        else:
            st.markdown("- US prevalence not in the curated map; for an orphan condition this is expected.")

        if row.get("combo_partner_1_name"):
            st.markdown("**Combination companion**")
            c1 = row["combo_partner_1_name"]
            syn1 = row.get("combo_partner_1_synergy")
            bridge1 = row.get("combo_partner_1_bridge_target")
            st.markdown(
                f"- Companion: **{c1}** (bridge target {bridge1}, synergy "
                f"{float(syn1):.2f}). The combination covers a disease target the primary substance misses."
            )
            if row.get("combo_partner_2_name"):
                st.markdown(
                    f"- Alternate: **{row['combo_partner_2_name']}** (bridge "
                    f"{row.get('combo_partner_2_bridge_target')}, "
                    f"synergy {float(row['combo_partner_2_synergy']):.2f})."
                )

        # Bundled PDF brief: inline iframe (sandboxed -- always renders) AND
        # download button. iframe lets the user read the brief in place;
        # download button is the fallback for browsers that block PDF iframes.
        substance_name = row.get("substance_name") or row.get("drug_name") or ""
        disease_name = row.get("disease_name") or ""
        brief_p = brief_path_for(substance_name, disease_name)
        if brief_p is not None:
            st.markdown("**One-page PDF brief**")
            pdf_b64 = base64.b64encode(brief_p.read_bytes()).decode("ascii")
            brief_html = (
                "<!doctype html><html><body style='margin:0'>"
                f"<iframe src='data:application/pdf;base64,{pdf_b64}' "
                "style='width:100%;height:680px;border:1px solid #e0e0e0;"
                "border-radius:6px' title='REPRISE brief'></iframe>"
                "</body></html>"
            )
            st.components.v1.html(brief_html, height=720, scrolling=False)
            st.download_button(
                label=f":page_facing_up: Download brief ({substance_name} -> {disease_name})",
                data=brief_p.read_bytes(),
                file_name=brief_p.name,
                mime="application/pdf",
                key=f"dl_{brief_p.name}",
            )
        else:
            st.caption("_No bundled PDF brief for this hypothesis -- regenerate via "
                       "`scripts/build_deal_memos.py --no-kol` to add one._")

        with st.expander("Raw row (all 36 columns)"):
            st.dataframe(row.to_frame().T, use_container_width=True, hide_index=True)

    # ------------------------------------------------------------------
    # Resolve the selected row from either the scatter OR the table.
    # Both widgets share the same head(2000) view. customdata in plotly
    # carries the dataframe index so we can resolve it directly.
    # ------------------------------------------------------------------
    selected_index = st.session_state.get("picker_index")
    if selected_index is None:
        sc_pts = (scatter_event.selection or {}).get("points") if hasattr(scatter_event, "selection") else []
        if sc_pts:
            p = sc_pts[0]
            pos = p.get("point_index") if "point_index" in p else p.get("pointIndex")
            if pos is not None and pos < len(scatter_view):
                selected_index = int(scatter_view.iloc[pos]["_df_idx"])
    if selected_index is None:
        sel = (event.selection or {}).get("rows") if hasattr(event, "selection") else []
        if sel:
            selected_index = view.head(2000).iloc[sel[0]].name

    if selected_index is not None:
        st.divider()
        st.markdown("## :pill: Hypothesis detail")
        render_brief(df.loc[selected_index])
    else:
        st.info("Pick from the dropdown above, click any scatter point, or click "
                "any table row to see the full mechanism / IP / market brief with "
                "chemical structure and inline PDF.")


# ======================================================================
# Tab 2 -- History & how REPRISE adds value
# ======================================================================
with tab_history:
    st.header("Drug repurposing: how the industry got here")
    st.markdown(
        "Most landmark repurposings of the past forty years were **serendipitous** "
        "-- a clinician noticed a side effect (sildenafil's vasodilation, "
        "propranolol's effect on a haemangioma, minoxidil's hair-growth side "
        "effect), an off-label experiment generated a positive trial signal "
        "(metformin in PCOS, naltrexone in alcohol dependence), or a pre-clinical "
        "screen revealed an unexpected cross-class activity (imatinib's KIT "
        "binding, bortezomib's mantle-cell-lymphoma response). The list below "
        "is a curated set of 54 historic wins drawn from FDA / EMA approvals and "
        "well-established off-label uses, used both to validate REPRISE and to "
        "make the historical pattern visible."
    )

    history = load_repurposing_history()
    if history.empty:
        st.info("Repurposing history dataset not bundled with this dashboard.")
    else:
        # ------------------------------------------------------------------
        # Timeline
        # ------------------------------------------------------------------
        st.subheader("Timeline of 54 historic repurposings")
        timeline = history.assign(
            decade=(history["year"] // 10 * 10).astype(int).astype(str) + "s",
            label=history["drug"] + " -> " + history["repurposed"].str.lower(),
        ).sort_values("year")
        fig_t = px.scatter(
            timeline,
            x="year",
            y=[1] * len(timeline),
            color="decade",
            hover_data={
                "drug": True, "original": True, "repurposed": True,
                "targets": True, "notes": True,
                "year": True, "decade": False,
            },
            labels={"x": "Year repurposed", "year": "Year repurposed"},
        )
        fig_t.update_traces(marker=dict(size=12, opacity=0.85,
                                         line=dict(width=0.6, color="#333")))
        fig_t.update_yaxes(visible=False, range=[0.8, 1.2])
        fig_t.update_layout(
            height=320,
            xaxis=dict(dtick=5, showgrid=True, gridcolor="#eaeaea"),
            legend=dict(orientation="h", yanchor="bottom", y=1.02,
                        xanchor="right", x=1, title=""),
            margin=dict(l=10, r=10, t=40, b=40),
            plot_bgcolor="white",
        )
        st.plotly_chart(fig_t, use_container_width=True)
        st.caption(
            "Hover any point for the original indication, the repurposed use, "
            "and the mechanism targets. The pace accelerates from the late 1990s "
            "onwards -- mirroring the rise of systematic target-disease databases "
            "(Open Targets, ChEMBL) that REPRISE now leverages."
        )

        # ------------------------------------------------------------------
        # Disease-disease network bridged by repurposed drugs
        # ------------------------------------------------------------------
        st.subheader("Disease network linked by repurposed drugs")
        st.markdown(
            "Each **node is an indication**, each **edge is a drug** that "
            "was repurposed from the original to the new use. Hub nodes "
            "(hypertension, epilepsy, rheumatoid arthritis, schizophrenia, "
            "depression) are the most fertile sources of repurposings -- a "
            "structural observation REPRISE captures and extends "
            "across the full disease ontology, not just the historic set."
        )

        import networkx as nx
        G = nx.Graph()
        for _, r in history.iterrows():
            o = r["original"].split(" / ")[0].split(",")[0].strip().lower()
            d = r["repurposed"].strip().lower()
            if not o or not d or o == d:
                continue
            G.add_edge(o, d, drug=r["drug"], year=int(r["year"]))

        if G.number_of_edges() == 0:
            st.info("Network has no usable edges from the bundled YAML.")
        else:
            pos = nx.spring_layout(G, k=1.1, seed=7, iterations=80)
            # edges
            edge_x, edge_y, edge_hover = [], [], []
            for u, v, attrs in G.edges(data=True):
                x0, y0 = pos[u]; x1, y1 = pos[v]
                edge_x += [x0, x1, None]
                edge_y += [y0, y1, None]
                edge_hover.append(f"{attrs['drug']} ({attrs['year']}): {u} -> {v}")
            edge_trace = go.Scatter(
                x=edge_x, y=edge_y, mode="lines",
                line=dict(width=1.0, color="#bbb"),
                hoverinfo="skip",
                showlegend=False,
            )
            # nodes -- size & color by degree (hubs stand out)
            node_x, node_y, node_text, node_hover, node_size, node_color = [], [], [], [], [], []
            for n in G.nodes():
                x, y = pos[n]
                node_x.append(x); node_y.append(y)
                deg = G.degree(n)
                node_size.append(10 + 6 * deg)
                node_color.append(deg)
                # label only the hubs (degree>=2) to avoid clutter
                node_text.append(n if deg >= 2 else "")
                drug_list = ", ".join(sorted({d["drug"] for _, _, d in G.edges(n, data=True)}))
                node_hover.append(
                    f"<b>{n}</b><br>degree (n drugs bridging it): {deg}"
                    f"<br>drugs: {drug_list}"
                )
            node_trace = go.Scatter(
                x=node_x, y=node_y, mode="markers+text",
                text=node_text, textposition="top center",
                textfont=dict(size=10),
                hovertext=node_hover, hoverinfo="text",
                marker=dict(
                    size=node_size,
                    color=node_color,
                    colorscale="Plasma",
                    showscale=True,
                    colorbar=dict(title="Drugs bridging<br>this indication",
                                  thickness=12, len=0.6),
                    line=dict(width=1, color="#222"),
                ),
                showlegend=False,
            )
            fig_n = go.Figure(data=[edge_trace, node_trace])
            fig_n.update_layout(
                height=560,
                showlegend=False,
                hovermode="closest",
                margin=dict(l=10, r=10, t=20, b=10),
                plot_bgcolor="white",
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
            )
            st.plotly_chart(fig_n, use_container_width=True)
            st.caption(
                "Hub indications (highlighted text, brighter colour) are sources "
                "from which many drugs have been repurposed. Hypertension alone "
                "spawned propranolol -> haemangioma & PTSD, eplerenone & "
                "spironolactone & carvedilol -> heart failure, "
                "verapamil -> migraine, prazosin -> PTSD, and minoxidil -> "
                "alopecia."
            )

    st.divider()
    st.subheader("How REPRISE is different and what it adds")
    cmp1, cmp2 = st.columns(2)
    with cmp1:
        st.markdown("**Historic repurposing**")
        st.markdown(
            "- **Serendipitous.** A clinician noticed a side effect, an "
            "off-label trial got lucky, a screen surfaced a cross-class hit.\n"
            "- **Slow.** Median **decades** between original use and "
            "FDA / EMA repurposed approval.\n"
            "- **One drug at a time.** Each story was its own narrative.\n"
            "- **Hit-driven only.** No systematic measurement of how often "
            "the mechanistic case actually shows up in the literature; no "
            "filter against crowded patent space; no orphan-prevalence "
            "overlay.\n"
            "- **No combination layer.** The BRAF + MEK doublet took years "
            "of separate oncology trials to identify."
        )
    with cmp2:
        st.markdown("**REPRISE adds**")
        st.markdown(
            "- **Systematic.** 3,996 approved drugs x 28,198 diseases "
            "= 176,272 ranked hypotheses, recomputed in ~12 minutes on a "
            "laptop.\n"
            "- **Measured precision.** Validated against this same 54-case "
            "set: **89% recovery (95% CI 77-96%)** at mech-support 0.30.\n"
            "- **Mechanism-aware.** Scoring is a noisy-OR over drug-target x "
            "target-disease evidence -- not a pattern match.\n"
            "- **Layered for action.** Directionality (does the drug push the "
            "target the way OT genetics predicts is therapeutic?), tissue, "
            "pathway, IP runway, orphan-prevalence, literature & patent prior.\n"
            "- **Combination companion finder.** Re-derives the BRAF + MEK "
            "doublet *de novo* in cardiofaciocutaneous syndrome -- the same "
            "combo already standard-of-care for BRAF V600E melanoma, without "
            "prior cross-reference.\n"
            "- **Operational output.** Per-hit one-page PDF brief that "
            "translates a row of numbers into a narrative suitable for "
            "clinical-collaborator outreach."
        )

    st.divider()
    st.markdown(
        "The point is **not to replace clinician judgement** but to compress "
        "what used to be a decades-long serendipity loop into a recomputable "
        "screen with measured precision and explicit regulatory / economic "
        "context. Use the Browse tab to interrogate the 176,272-row output; "
        "the FAQ tab for what every number means."
    )


# ======================================================================
# Tab 3 -- FAQ
# ======================================================================
with tab_faq:
    st.header("FAQ -- how to read this dashboard")
    st.caption("Skim the headline first. Then use the column-by-column reference "
               "to decode any number you see in the Browse tab.")

    # --- TLDR -----------------------------------------------------------
    st.markdown("### What is REPRISE in one paragraph?")
    st.markdown(
        "REPRISE is a screen that asks, for every approved drug and every disease "
        "in the public ontology: *if a clinical collaborator wanted to repurpose "
        "this drug for this disease, how strong is the mechanistic case, how novel "
        "is the indication, how is the IP runway, how large is the patient "
        "population, and how crowded is the literature?* It produces a single "
        "**opportunity score** per (drug, disease) pair across **176,272 hypotheses**, "
        "drawing on Open Targets, ChEMBL, Reactome, STRING, EFO, FDA Orange Book, "
        "Orphanet, PubMed, Europe PMC, ClinicalTrials.gov, and Lens.org. It was "
        "validated by recovering **48 of 54 known repurposing wins (89%, 95% CI "
        "77-96%)** at a mechanism-support threshold of 0.30."
    )

    st.markdown("### What's the main message of the paper?")
    st.markdown(
        "Mechanism-only ranking is necessary but not sufficient -- the layers that "
        "separate an actionable hypothesis from a plausible one are **directionality** "
        "(does the drug push the target the way the disease's genetics says is "
        "therapeutic?), **tissue expression**, **pathway context**, **IP runway**, "
        "**market size**, and **literature crowding**. The paper shows the screen "
        "rediscovers the BRAF + MEK doublet *de novo* in cardiofaciocutaneous "
        "syndrome (the same combination already standard-of-care in BRAF V600E "
        "metastatic melanoma) and recapitulates metformin's mitochondrial complex-I "
        "mechanism for polycystic ovary syndrome, not the AMPK route the literature "
        "most often cites."
    )

    st.markdown("### How do I use the dashboard?")
    st.markdown(
        "- **Browse tab** -- start with the sidebar filters. Type a drug or disease "
        "name, drag the Opportunity slider up, optionally tick *Orphan indications* "
        "or *Has Lens patent activity*.\n"
        "- **Click any row** -- the bottom panel renders the full per-hypothesis "
        "brief: mechanism rationale, IP / literature / patent context, clinical "
        "opportunity with prevalence source attribution, and a combination "
        "companion if one was found.\n"
        "- **Hover the column headers** in the table for inline definitions. The "
        "definitions of every score are also expanded below."
    )

    st.divider()
    st.markdown("### How are the numbers computed?")
    st.caption("One expander per column. Click to read the definition.")

    with st.expander("Active ingredient -- what's in this column?"):
        st.markdown(
            "ChEMBL ships drugs at both the parent-molecule level (e.g. SILDENAFIL, "
            "`CHEMBL192`) and the salt / formulation level (SILDENAFIL CITRATE, "
            "`CHEMBL1737`). Mechanism-of-action curation can live on either level. "
            "REPRISE rolls all formulations up to the **`molecule_hierarchy` parent** "
            "and unions targets across them, so each row in this column is the "
            "canonical active ingredient. The detail panel's *Raw row* expander "
            "shows the underlying ChEMBL ID."
        )

    with st.expander("Mech (mechanism support) -- the core score"):
        st.markdown(
            "$\\mathrm{mech\\_support} = 1 - \\prod_i (1 - w_i \\times \\mathrm{assoc}_i)$\n\n"
            "A **noisy-OR** over the drug's direct targets (weight 1.0) and "
            "optionally STRING-neighbour targets (weight 0.5 x STRING confidence). "
            "For each (drug_target, disease) edge the contribution is "
            "`target_weight x OT_association_score`. Open Targets associations "
            "below 0.1 are dropped to suppress noise. Value range **[0, 1]**: 0 = no "
            "mechanism overlap; 1 = saturated."
        )

    with st.expander("Novelty"):
        st.markdown(
            "1 if the (drug, disease) pair is **not** on ChEMBL's known indication "
            "list for that drug, and **not** present for any other drug sharing the "
            "exact same set of (direct target, action_type) pairs (target-class "
            "rollup). The screen walks one EFO ontology hop outward to catch "
            "close-but-distinct conditions. STRING neighbours are excluded from the "
            "class keys so weak-edge targets don't inflate the rollup."
        )

    with st.expander("Direction (direction_status / direction_factor)"):
        st.markdown(
            "From Open Targets genetic evidence, restricted to seven curated "
            "human-genetics sources (`ot_genetics_portal`, `gene_burden`, `eva`, "
            "`gene2phenotype`, `clingen`, `genomics_england`, `orphanet`). IMPC "
            "mouse-knockout evidence is excluded because LoF in the mouse rarely "
            "corresponds to therapeutic direction in human monogenic disease.\n\n"
            "Combining `variantEffect in {LoF, GoF}` with `directionOnTrait in "
            "{risk, protective}` gives a therapeutic direction in {+1, 0, -1} per "
            "(target, disease). If the drug's action class is compatible with the "
            "predicted direction (e.g. an antagonist for a GoF/risk target), the "
            "direction factor multiplies opportunity by 1; if opposed, by 0.5; if "
            "unknown, by 1.0 (no penalty)."
        )

    with st.expander("Tissue (tissue_status / tissue_factor)"):
        st.markdown(
            "Per (target, disease), score is the max over the disease's relevant "
            "tissues of (`disease_relevance x target_expression`), bucketed:\n"
            "- **expressed** (>= 0.25): tissue factor 1.0\n"
            "- **low**: linear ramp 0.3 -> 1.0\n"
            "- **absent** (measured zero): tissue factor 0.3\n"
            "- **unknown** (no measurement): tissue factor 0.7\n\n"
            "Multiplicative on opportunity. The 0.7-not-1.0 default for *unknown* "
            "is deliberate -- it modestly penalises hypotheses where we genuinely "
            "have no expression data, without erasing them."
        )

    with st.expander("Phylo (phylo_score / phylo_factor)"):
        st.markdown(
            "From Open Targets evidence of type `animal_model` (default source: "
            "**IMPC** via PhenoDigm), the max `phylo_score` across the drug's "
            "targets per (drug, disease). `phylo_factor = 1 + 0.5 x phylo_score`. "
            "Asymmetric: presence boosts opportunity, absence does not penalise."
        )

    with st.expander("Pathway (pathway_score / n_pathway_overlap)"):
        st.markdown(
            "Reactome co-membership between the drug's targets and the disease's "
            "strongly associated targets, restricted to pathways with <= 80 genes "
            "(specificity filter). Only **indirect** bridges count -- the drug's "
            "target must differ from the disease target -- so direct-target hits "
            "are not double-counted. "
            "$\\mathrm{pathway\\_factor} = 1 + 0.3 \\times \\log_1\\!\\!p(n\\_overlap) "
            "/ \\log_1\\!\\!p(5)$, capped at 5 overlaps."
        )

    with st.expander("Opportunity -- the composite score"):
        st.markdown(
            "$\\mathrm{opportunity} = \\mathrm{mech} \\times \\mathrm{novelty} \\times "
            "\\mathrm{direction} \\times \\mathrm{tissue} \\times \\mathrm{phylo} \\times "
            "\\mathrm{pathway} \\,/\\, \\sqrt{n\\_drug\\_targets}$\n\n"
            "The `sqrt(n_targets)` denominator is a mild penalty against highly "
            "promiscuous drugs whose mech_support might be inflated by spurious "
            "weak edges. After this score, the substance-grouping step rolls "
            "formulations up to the active ingredient, and the literature pass "
            "damps by `(1 - 0.5 x investigation_prior)` for the top 5,000 rows."
        )

    with st.expander("Investigation prior (inv_prior) and the patent column"):
        st.markdown(
            "For the top-N hypotheses (default 5,000), REPRISE queries PubMed, "
            "Europe PMC, ClinicalTrials.gov and Lens.org for the (target x disease) "
            "combination. Each source's count is log-saturated against a per-source "
            "saturation point (e.g. 200 PubMed papers saturates to 1.0). The four "
            "log-counts are weighted-averaged into `investigation_prior in [0, 1]`. "
            "Opportunity is damped by `(1 - 0.5 x investigation_prior)` -- a "
            "well-investigated hypothesis falls in rank.\n\n"
            "Lens.org coverage can be extended past the top-5,000 cap by setting "
            "`literature.lens_top_n` in the config -- useful for catching crowded "
            "patent frontiers that the cheap sources don't surface."
        )

    with st.expander("Severity damping (severity_concern column)"):
        st.markdown(
            "Mechanism scoring has one recurring failure mode in monogenic disease: "
            "if a receptor is **broken** (loss of function) and the drug is an "
            "**agonist**, the drug cannot rescue the broken protein. The screen "
            "flags such hypotheses when the disease name carries severity language "
            "(`deficiency`, `complete absence`, `Donohue syndrome`, etc.) AND a "
            "gene named in the disease is also a direct target of the drug AND the "
            "drug's action is agonistic. Flagged rows are tagged "
            "`severe_loF_agonist` and their opportunity is reduced by 70%. The "
            "*Hide severity-flagged* sidebar checkbox suppresses them by default. "
            "Adjacent partial-LoF rescuable cases (TZD -> PPARG lipodystrophy) "
            "are preserved."
        )

    with st.expander("IP runway columns (Orange Book)"):
        st.markdown(
            "From the FDA Orange Book (Approved Drug Products with Therapeutic "
            "Equivalence Evaluations), parsed per active ingredient:\n"
            "- **Patent yr (`latest_patent_year`)** -- latest listed patent expiry\n"
            "- **`loe_year`** -- max of patent and exclusivity expiry; this is the "
            "year the drug becomes fair game for generic competition\n"
            "- **Generic (`has_generic`)** -- True if any ANDA generic has been "
            "filed against the active ingredient\n\n"
            "The Orange Book is small-molecule only. Biologics (mAbs, peptides) "
            "live in the Purple Book and currently show empty Orange Book signals "
            "in this screen -- the per-hit brief surfaces this as `check the "
            "Purple Book`."
        )

    with st.expander("Combination partner & synergy"):
        st.markdown(
            "For each top-N primary hypothesis, the companion finder evaluates "
            "candidate combination substances whose direct targets bridge **strong "
            "disease targets the primary does not hit**. Synergy is "
            "`combo_mech_support - primary_mech_support` under the noisy-OR. A "
            "target-overlap filter rejects same-class redundancy (>50% target "
            "overlap with the primary is rejected). The Combo partner column "
            "shows the highest-synergy companion; click the row to see the bridge "
            "target and any alternate companion."
        )

    with st.expander("Orphan flag and US patients"):
        st.markdown(
            "`is_orphan = True` when the curated US patient count is below the "
            "**200,000** threshold from the FDA Orphan Drug Act, which qualifies "
            "the indication for 7-year orphan exclusivity and premium pricing. "
            "Patient counts come from two sources: a 122-disease curated CSV "
            "drawn from CDC, SEER, and disease-foundation estimates, plus Orphanet "
            "structured prevalence for ~5,100 rare diseases bridged to OT MONDO "
            "via dbXRefs."
        )

    st.divider()
    st.markdown("### A few honest caveats")
    st.markdown(
        "- The screen is a **research artefact**. Scores are computational "
        "priors, not clinical recommendations.\n"
        "- **Open Targets coverage is the precision ceiling.** The six backtest "
        "misses (minoxidil/alopecia, propranolol/hemangioma, acetazolamide/IIH, "
        "verapamil/migraine, anakinra/CAPS, eculizumab/aHUS) all trace to OT "
        "associations missing the canonical mechanism target. They are not "
        "pipeline-logic errors and will lift automatically with each OT release.\n"
        "- The literature pass currently covers the **top 5,000** hypotheses by "
        "default; rows beyond that have NA in the literature columns.\n"
        "- The screen does not run drug-drug interaction or safety checks on "
        "combination companions -- treat the synergy column as a hypothesis "
        "generator, not a recommendation."
    )

    st.divider()
    st.markdown(
        f"Full manuscript and source code: "
        f"[github.ugent.be/wvcrieki/repurposed]({PAPER_URL})"
    )


st.divider()
st.caption(
    "REPRISE is a research artefact, not a clinical recommendation. "
    "All scores are computational priors. "
    f"[Source]({PAPER_URL})  |  manuscript: `manuscript/repurposing_screen_manuscript.pdf`"
)
