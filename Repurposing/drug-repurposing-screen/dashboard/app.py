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

from pathlib import Path

import pandas as pd
import streamlit as st

DATA = Path(__file__).parent / "data" / "repurposing_hypotheses.parquet"
PAPER_URL = "https://github.ugent.be/wvcrieki/repurposed"

st.set_page_config(
    page_title="REPRISE -- drug-repurposing screen",
    page_icon=":pill:",
    layout="wide",
)


# ----------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def load_data() -> pd.DataFrame:
    df = pd.read_parquet(DATA)
    # Coerce a few columns to friendly types for filtering
    for col in ("opportunity", "mechanistic_support", "novelty",
                "investigation_prior", "pathway_score"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
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
c5.metric("Backtest precision", "89%", help="48/54 known repurposings recovered, 95% CI 77-96%")


# ----------------------------------------------------------------------
# Sidebar filters
# ----------------------------------------------------------------------
st.sidebar.header("Filters")

drug_q = st.sidebar.text_input("Drug name contains").strip().upper()
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


# ----------------------------------------------------------------------
# Filter application
# ----------------------------------------------------------------------
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


# ----------------------------------------------------------------------
# Results table
# ----------------------------------------------------------------------
st.subheader(f"Ranked hypotheses ({len(view):,} match{'es' if len(view) != 1 else ''})")

if view.empty:
    st.info("No hypotheses match these filters. Loosen them in the sidebar.")
    st.stop()

display_cols = [
    c for c in [
        "rank", "substance_name", "disease_name", "lead_target",
        "opportunity", "mechanistic_support", "novelty", "direction_status",
        "tissue_status", "is_orphan", "us_patients",
        "latest_patent_year", "has_generic",
        "pubmed_count", "trial_count", "patent_count", "investigation_prior",
        "combo_partner_1_name", "combo_partner_1_synergy",
    ] if c in view.columns
]

# st.dataframe with row selection (Streamlit 1.31+ pattern)
event = st.dataframe(
    view[display_cols].head(2000),
    use_container_width=True,
    hide_index=True,
    column_config={
        "rank": st.column_config.NumberColumn("rank", format="%d"),
        "opportunity": st.column_config.NumberColumn(format="%.3f"),
        "mechanistic_support": st.column_config.NumberColumn("mech", format="%.3f"),
        "novelty": st.column_config.NumberColumn(format="%.2f"),
        "investigation_prior": st.column_config.NumberColumn("inv_prior", format="%.2f"),
        "combo_partner_1_synergy": st.column_config.NumberColumn("combo_syn", format="%.3f"),
        "us_patients": st.column_config.NumberColumn("US patients", format="%d"),
    },
    on_select="rerun",
    selection_mode="single-row",
    key="hits_table",
)

st.caption("Showing the top 2,000 rows after filtering. Refine in the sidebar to "
           "drill in. Click any row for the full per-hypothesis brief.")


# ----------------------------------------------------------------------
# Per-hit detail panel -- the same content as the PDF brief
# ----------------------------------------------------------------------
def _fmt_int(v) -> str:
    try:
        return f"{int(float(v)):,}"
    except Exception:
        return "n/a"


def render_brief(row: pd.Series) -> None:
    substance = str(row.get("substance_name") or row.get("drug_name") or "?")
    disease = str(row.get("disease_name") or "?")
    st.subheader(f"{substance} -> {disease}")

    # Why this hypothesis is real
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

    # IP runway + literature / patent prior
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

    # Clinical opportunity
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

    # Combination companion
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

    with st.expander("Raw row (all 36 columns)"):
        st.dataframe(row.to_frame().T, use_container_width=True, hide_index=True)


sel = (event.selection or {}).get("rows") if hasattr(event, "selection") else []
if sel:
    st.divider()
    selected_index = view.head(2000).iloc[sel[0]].name
    render_brief(df.loc[selected_index])
else:
    st.info("Click any row above to see the full mechanism / IP / market brief.")


st.divider()
st.caption(
    "REPRISE is a research artefact, not a clinical recommendation. "
    "All scores are computational priors. "
    f"[Source]({PAPER_URL})  |  manuscript: `manuscript/repurposing_screen_manuscript.pdf`"
)
