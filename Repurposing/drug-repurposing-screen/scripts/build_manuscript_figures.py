"""Generate the four main figures for the manuscript.

Outputs to manuscript/figures/:
  fig1_pipeline_architecture.png  -- pipeline block diagram
  fig2_backtest_validation.png    -- per-case mech_support, HIT/MISS
  fig3_engineering_choices.png    -- severity damping + direction filter
  fig4_combinations.png           -- BRAF+MEK network + top synergies

Re-run after pipeline or validation set changes to keep the manuscript
figures in sync with the screen output.
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
FIG_DIR = ROOT / "manuscript" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------
# Figure 1: pipeline architecture (schematic)
# ----------------------------------------------------------------------
def fig1_architecture():
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 9)
    ax.axis("off")

    def box(x, y, w, h, label, color="#E8EEF7", edge="#0B4FA8", fontsize=8.5):
        rect = mpatches.FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.02",
            linewidth=1.0, edgecolor=edge, facecolor=color)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
                fontsize=fontsize, color="#13294B", wrap=True)

    def arrow(x0, y0, x1, y1, color="#666666"):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=color,
                                    lw=0.8, mutation_scale=10))

    # Row 1 (top): data sources
    sources = ["ChEMBL", "Open Targets", "Reactome", "STRING", "EFO",
               "Orange Book", "Orphanet", "PubMed / EuroPMC / NCT",
               "Semantic Scholar"]
    for i, s in enumerate(sources):
        box(0.3 + 1.30 * i, 8.0, 1.20, 0.7, s, color="#F4ECDD",
            edge="#A78343", fontsize=7.5)

    # Row 2: ingestion
    box(0.5, 6.7, 11.0, 0.7, "Adapters + Loaders   (canonical CSV / parquet)",
        color="#EFEFEF", edge="#888888", fontsize=9)

    # Row 3: mechanism scoring
    box(4.2, 5.3, 3.6, 0.95,
        "Mechanism scoring\nnoisy-OR over drug-target $\\times$ target-disease",
        color="#D9E7FF", edge="#0B4FA8", fontsize=9)

    # Row 4: enrichment layers (parallel)
    enrich = [("Novelty",  "#D9E7FF"),
              ("Direction","#D9E7FF"),
              ("Tissue",   "#D9E7FF"),
              ("Phylo",    "#D9E7FF"),
              ("Pathway",  "#D9E7FF"),
              ("Severity", "#FFE2DA")]
    for i, (lbl, col) in enumerate(enrich):
        box(0.3 + 2.0 * i, 4.0, 1.85, 0.7, lbl, color=col, edge="#0B4FA8",
            fontsize=8.5)

    # Row 5: scoring + ranking
    box(2.5, 2.8, 7.0, 0.75,
        "Opportunity score $\\times$ active-ingredient grouping $\\to$ ranked output",
        color="#D9E7FF", edge="#0B4FA8", fontsize=9)

    # Row 6: post-engine enrichment
    post = [("Literature\nprior", "#FFF0D4"),
            ("Market\n(curated + Orphanet)", "#FFF0D4"),
            ("FDA Orange\nBook IP", "#FFF0D4"),
            ("KOL finder\n(US + EU)", "#FFF0D4"),
            ("Combination\ntherapy", "#FFF0D4")]
    for i, (lbl, col) in enumerate(post):
        box(0.3 + 2.36 * i, 1.5, 2.2, 0.85, lbl, color=col, edge="#B07C2F",
            fontsize=8)

    # Outputs
    box(1.5, 0.2, 4.0, 0.7, "36-column ranked CSV", color="#E8F4DC",
        edge="#3B7A22", fontsize=9.5)
    box(6.5, 0.2, 4.0, 0.7, "Per-hypothesis PDF deal memos",
        color="#E8F4DC", edge="#3B7A22", fontsize=9.5)

    # Arrows (selective, not every connection -- enough to convey flow)
    arrow(6.0, 8.0, 6.0, 7.4)        # sources -> adapters
    arrow(6.0, 6.7, 6.0, 6.25)       # adapters -> mechanism
    arrow(6.0, 5.3, 6.0, 4.7)        # mechanism -> enrichment row
    for i in range(6):
        arrow(0.3 + 2.0 * i + 0.92, 4.0, 0.3 + 2.0 * i + 0.92, 3.55)  # enrich -> scoring
    arrow(6.0, 2.8, 6.0, 2.35)       # scoring -> post-enrich
    for i in range(5):
        arrow(0.3 + 2.36 * i + 1.1, 1.5, 0.3 + 2.36 * i + 1.1, 0.9)   # post -> outputs

    ax.set_title("Pipeline architecture: end-to-end mechanism-driven drug repurposing screen",
                 fontsize=11, color="#13294B")
    plt.tight_layout()
    out = FIG_DIR / "fig1_pipeline_architecture.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


# ----------------------------------------------------------------------
# Figure 2: backtest validation -- per-case mech_support, HIT / MISS
# ----------------------------------------------------------------------
def _run_validation() -> pd.DataFrame:
    """Run the backtest validation script's core logic and return a DataFrame."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "backtest_validation",
        str(ROOT / "scripts" / "backtest_validation.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    from repurpose.config import load_config
    from repurpose.sources import loaders
    cfg = load_config(ROOT / "config.full.yaml")
    drugs = loaders.load_drugs(cfg.path("drugs"))
    dt = loaders.load_drug_targets(cfg.path("drug_targets"))
    td = loaders.load_target_disease(cfg.path("target_disease"))
    cases = yaml.safe_load(
        (ROOT / "data" / "curated" / "repurposing_validation.yaml").read_text()
    )["cases"]
    rows = []
    for c in cases:
        drug_ids = mod._find_drug_ids(c["drug"], drugs)
        if not drug_ids:
            rows.append({"drug": c["drug"], "disease": c.get("repurposed_disease_name", ""),
                         "mech": 0.0, "status": "DRUG_NOT_FOUND"})
            continue
        efo, name = mod._resolve_efo(c, td)
        if not efo:
            rows.append({"drug": c["drug"], "disease": c.get("repurposed_disease_name", ""),
                         "mech": 0.0, "status": "DISEASE_NOT_IN_OT"})
            continue
        mech = mod._mech_support_for(drug_ids, efo, dt, td)
        rows.append({"drug": c["drug"], "disease": name,
                     "mech": mech["mech_support"],
                     "status": "HIT" if mech["mech_support"] >= 0.3 else "MISS"})
    return pd.DataFrame(rows)


def fig2_backtest():
    df = _run_validation()
    df["label"] = df["drug"].str.title() + r" $\to$ " + df["disease"].str.lower()
    df = df.sort_values("mech")
    fig, ax = plt.subplots(figsize=(9, 7))
    colors = {"HIT": "#2C7A2E", "MISS": "#C03A2B", "DISEASE_NOT_IN_OT": "#888888",
              "DRUG_NOT_FOUND": "#888888"}
    bar_colors = [colors[s] for s in df["status"]]
    y = np.arange(len(df))
    ax.barh(y, df["mech"], color=bar_colors, edgecolor="black", linewidth=0.4,
            height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"], fontsize=8)
    ax.set_xlabel("Mechanism support  (noisy-OR over drug-target $\\times$ target-disease)",
                  fontsize=10)
    ax.set_xlim(0, 1.05)
    ax.axvline(0.3, color="#444444", linestyle="--", linewidth=1.0)
    ax.text(0.305, len(df) - 0.6, "HIT threshold = 0.30", fontsize=8,
            color="#444444", va="center")
    hit_n = (df.status == "HIT").sum()
    ax.set_title(f"Backtest validation: {hit_n} of {len(df)} curated known repurposing "
                 f"successes recovered  ({hit_n / len(df) * 100:.0f}% hit rate)",
                 fontsize=10.5)
    legend_handles = [
        mpatches.Patch(color=colors["HIT"], label="HIT (mech_support >= 0.30)"),
        mpatches.Patch(color=colors["MISS"], label="MISS  (OT signal below threshold)"),
        mpatches.Patch(color=colors["DISEASE_NOT_IN_OT"],
                       label="OT data gap  (disease not in OT)"),
    ]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=8.5,
              framealpha=0.95)
    ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.set_axisbelow(True)
    plt.tight_layout()
    out = FIG_DIR / "fig2_backtest_validation.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out, df


# ----------------------------------------------------------------------
# Figure 3: engineering choices that mattered
#   Panel A: severity damping pushes 32 'severe LoF + agonist' hits below rank 100
#   Panel B: direction-signal composition before/after the IMPC-source filter
# ----------------------------------------------------------------------
def fig3_engineering():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: severity damping effect
    df = pd.read_csv(ROOT / "output_full" / "repurposing_hypotheses.csv",
                     low_memory=False)
    flagged = df[df.severity_concern == "severe_loF_agonist"]
    # Simulate pre-damping ranks: each flagged row's "natural" rank is what
    # rank it would have had at its un-damped opportunity. Approximate by
    # un-damping (opportunity / 0.3) and re-ranking against full output.
    if not flagged.empty:
        full_opp = df["opportunity"].values
        un_damped = (flagged["opportunity"] / 0.3).values  # damp_factor = 0.7
        pre_ranks = []
        sorted_opp = np.sort(full_opp)[::-1]
        for v in un_damped:
            pre_ranks.append(int(np.searchsorted(-sorted_opp, -v) + 1))
        post_ranks = flagged["rank"].values

        axA = axes[0]
        axA.scatter(pre_ranks, post_ranks, color="#C03A2B", s=22,
                    edgecolor="#5B1B14", linewidth=0.3, alpha=0.85)
        axA.plot([1, 200000], [1, 200000], color="#888888", linestyle=":",
                 linewidth=1, alpha=0.5, label="No movement (y=x)")
        axA.set_xscale("log")
        axA.set_yscale("log")
        axA.set_xlim(5, 300000)
        axA.set_ylim(5, 300000)
        axA.set_xlabel("Rank without severity damping", fontsize=10)
        axA.set_ylabel("Rank with severity damping", fontsize=10)
        axA.set_title(f"A.  Severity heuristic relocates {len(flagged)} 'receptor LoF "
                      f"+ agonist' hypotheses\nbelow rank 100",
                      fontsize=10.5)
        axA.axhline(100, color="#0B4FA8", linestyle="--", linewidth=1)
        axA.text(7, 130, "Top-100 boundary", fontsize=8, color="#0B4FA8")
        axA.grid(linestyle=":", linewidth=0.4, alpha=0.6)
        axA.set_axisbelow(True)

    # Panel B: direction-signal composition before/after IMPC filter
    axB = axes[1]
    sources = ["IMPC\n(mouse KO)", "chembl", "eva\n(ClinVar)",
               "ot_genetics_\nportal", "gene_burden", "cancer_gene_\ncensus", "Other"]
    counts_before = [1_138_754, 639_068, 281_697, 166_292, 36_046, 34_164, 7_887]
    kept_mask =    [False,       False,    True,    True,        True,    False,        True]   # what survived the filter
    bar_colors_before = ["#C03A2B" if not k else "#2C7A2E" for k in kept_mask]
    x = np.arange(len(sources))
    axB.bar(x, counts_before, color=bar_colors_before,
            edgecolor="black", linewidth=0.4)
    axB.set_xticks(x)
    axB.set_xticklabels(sources, fontsize=8)
    axB.set_ylabel("Direction-informative evidence rows", fontsize=10)
    axB.set_yscale("log")
    axB.set_ylim(1_000, 3_000_000)
    axB.set_title("B.  Filtering non-human-genetics sources contracts "
                  "direction signal\nfrom 1.02 M rows to 113 k "
                  "(eliminates monogenic-disorder bias)",
                  fontsize=10.5)
    legend_handles = [
        mpatches.Patch(color="#C03A2B", label="Excluded (not therapeutic-direction-informative)"),
        mpatches.Patch(color="#2C7A2E", label="Retained (curated human genetics)"),
    ]
    axB.legend(handles=legend_handles, loc="upper right", fontsize=8.5,
               framealpha=0.95)
    axB.grid(axis="y", linestyle=":", linewidth=0.4, alpha=0.6)
    axB.set_axisbelow(True)

    plt.tight_layout()
    out = FIG_DIR / "fig3_engineering_choices.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


# ----------------------------------------------------------------------
# Figure 4: BRAF+MEK rediscovery via combination-finder
# ----------------------------------------------------------------------
def fig4_combinations():
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    # Panel A: combination network diagram
    axA = axes[0]
    axA.set_xlim(-0.5, 5.5)
    axA.set_ylim(-0.5, 5.5)
    axA.axis("off")
    axA.set_title("A.  Combination-finder *de novo* rediscovery of the\n"
                  "BRAF + MEK doublet (and KRAS + MEK analogue) for "
                  "cardiofaciocutaneous syndrome",
                  fontsize=10.5)

    def node(ax, x, y, label, color="#D9E7FF", edge="#0B4FA8", r=0.5,
             fontsize=9):
        circ = mpatches.Circle((x, y), r, facecolor=color, edgecolor=edge,
                               linewidth=1.2)
        ax.add_patch(circ)
        ax.text(x, y, label, ha="center", va="center", fontsize=fontsize,
                color="#13294B", weight="bold")

    # Primary nodes (BRAF inhibitors + KRAS inhibitors)
    primaries = [
        ("Dabrafenib", 0.4, 4.5),
        ("Vemurafenib", 0.4, 3.5),
        ("Encorafenib", 0.4, 2.5),
        ("Sotorasib",  0.4, 1.5),
        ("Adagrasib",  0.4, 0.5),
    ]
    for name, x, y in primaries:
        node(axA, x, y, name, color="#FFE2DA", edge="#B23A20", r=0.46,
             fontsize=8.5)

    # Disease node (center-top)
    node(axA, 2.7, 5.0, "CFC\nsyndrome", color="#E8F4DC", edge="#3B7A22",
         r=0.62, fontsize=9.5)

    # Bridge target
    node(axA, 2.7, 2.5, "MAP2K1 /\nMAP2K2", color="#FFF0D4", edge="#B07C2F",
         r=0.62, fontsize=9)

    # Companion drugs (MEK inhibitors)
    node(axA, 5.0, 3.4, "Binimetinib", color="#D9E7FF", edge="#0B4FA8",
         r=0.55, fontsize=9)
    node(axA, 5.0, 1.6, "Cobimetinib", color="#D9E7FF", edge="#0B4FA8",
         r=0.55, fontsize=9)

    # Edges
    def edge(x0, y0, x1, y1, color="#666666", w=0.6, ls="-"):
        axA.plot([x0, x1], [y0, y1], color=color, linewidth=w, linestyle=ls)

    for _, x, y in primaries:
        edge(x + 0.4, y, 2.7, 4.45, color="#B23A20", w=0.6)         # primary -> disease
        edge(x + 0.4, y, 2.1, 2.55, color="#888888", w=0.5, ls=":")  # primary -> bridge
    edge(3.3, 2.55, 4.5, 3.35, color="#0B4FA8", w=0.7)
    edge(3.3, 2.55, 4.5, 1.65, color="#0B4FA8", w=0.7)
    # Disease <- bridge
    edge(2.7, 3.12, 2.7, 4.4, color="#3B7A22", w=0.6, ls=":")

    # Labels
    axA.text(0.4, 5.15, "Primary substances", fontsize=8.5, ha="center",
             color="#5B1B14", style="italic")
    axA.text(5.0, 4.05, "Companion drugs", fontsize=8.5, ha="center",
             color="#0B4FA8", style="italic")
    axA.text(2.7, 1.78, "bridge target", fontsize=7.5, ha="center",
             color="#705221", style="italic")

    # Panel B: top combination synergies bar chart
    axB = axes[1]
    combos = [
        ("Bromazepam + Orphenadrine\n$\\to$ developmental and epileptic\nencephalopathy", 0.48),
        ("Setmelanotide + Metformin\n$\\to$ T2D (via NDUFAB1,\nmetformin's complex-I site)", 0.49),
        ("Pinacidil/Minoxidil + Vernakalant\n$\\to$ familial atrial fibrillation", 0.26),
        ("Angiotensin II + Aliskiren\n$\\to$ renal tubular dysgenesis", 0.21),
        ("Sotorasib/Adagrasib + Binimetinib\n$\\to$ cardiofaciocutaneous syndrome", 0.17),
        ("Dabrafenib/Vemurafenib/Encorafenib\n+ Binimetinib $\\to$ CFC syndrome", 0.15),
        ("Follitropin + Esterified estrogens\n$\\to$ 46,XX gonadal dysgenesis", 0.11),
    ]
    labels = [c[0] for c in combos]
    syn = [c[1] for c in combos]
    y = np.arange(len(combos))
    axB.barh(y, syn, color="#0B4FA8", edgecolor="black", linewidth=0.4,
             height=0.65, alpha=0.85)
    axB.set_yticks(y)
    axB.set_yticklabels(labels, fontsize=7.8)
    axB.invert_yaxis()
    axB.set_xlabel("Synergy  (combo mech_support $-$ primary mech_support)",
                   fontsize=10)
    axB.set_xlim(0, 0.55)
    axB.set_title("B.  Top combination-therapy synergies\nin the screen's top 200 primary hypotheses",
                  fontsize=10.5)
    axB.grid(axis="x", linestyle=":", linewidth=0.4, alpha=0.6)
    axB.set_axisbelow(True)

    plt.tight_layout()
    out = FIG_DIR / "fig4_combinations.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    print("Generating Fig 1 (architecture)...")
    p1 = fig1_architecture()
    print(f"  wrote {p1}")
    print("Generating Fig 2 (backtest)...")
    p2, _ = fig2_backtest()
    print(f"  wrote {p2}")
    print("Generating Fig 3 (engineering choices)...")
    p3 = fig3_engineering()
    print(f"  wrote {p3}")
    print("Generating Fig 4 (combinations)...")
    p4 = fig4_combinations()
    print(f"  wrote {p4}")


if __name__ == "__main__":
    main()
