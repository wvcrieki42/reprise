"""Generate one-page PDF deal memos for KOL outreach on screen hits.

Each memo is a stand-alone PDF aimed at convincing a Key Opinion Leader
to start a three-way collaboration: the KOL leads the clinical /
translational side, our bioinformatics group brings analytics support,
and the substance's IP holder (visible via the Orange Book NDA holder)
supplies the compound under an investigator-initiated trial framework.

Narrative shape per memo:

    HYPOTHESIS         Drug -> Disease, one line.
    WHY THIS IS REAL   Mechanism (mech_support, novelty, direction, tissue,
                       phylogenetics, pathway), severity flag absence.
    WHY NOW            Orange Book IP runway, generic status, NDA holder.
                       Literature volume (open ground vs trodden).
    CLINICAL TARGET    US patients, orphan status, prevalence source.
    WHY YOU            KOL name, h-index, institution, query-specific pub
                       count; email if available.
    WHAT WE PROPOSE    Our group's bioinformatics pitch (from profile YAML),
                       suggested collaboration model.
    NEXT STEP          Concrete ask + signature block.

Defaults pull from data/curated/memo_profile.yaml. Override with
--profile, --pi-name, --pi-email, --group-name on the command line.

Usage:
    python scripts/build_deal_memos.py
        [--input output_full/repurposing_hypotheses.csv]
        [--output output_full/memos/]
        [--top 5]
        [--region us|eu|both]    # default: both, one memo per region with a KOL
        [--orphan-only]
        [--filter-min-opportunity 0.5]
"""
from __future__ import annotations
import argparse
import re
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PROFILE = ROOT / "data" / "curated" / "memo_profile.yaml"
DEFAULT_INPUT = ROOT / "output_full" / "repurposing_hypotheses.csv"
DEFAULT_OUTPUT = ROOT / "output_full" / "memos"


# ----------------------------------------------------------------------
# Profile loading
# ----------------------------------------------------------------------
def _load_profile(path: Path, overrides: dict) -> dict:
    profile = yaml.safe_load(Path(path).read_text())
    for k, v in overrides.items():
        if v is None:
            continue
        if k in ("pi_name", "pi_email", "group_name"):
            section, key = ("group",
                            {"pi_name": "pi_name", "pi_email": "pi_email",
                             "group_name": "name"}[k])
            profile.setdefault(section, {})[key] = v
    return profile


# ----------------------------------------------------------------------
# Narrative helpers
# ----------------------------------------------------------------------
def _fmt_int(x) -> str:
    if pd.isna(x):
        return "n/a"
    return f"{int(x):,}"


def _s(v) -> str:
    """Coerce a possibly-NaN-or-None value to a clean string."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    return str(v).strip()


def _fmt_pct(x: float) -> str:
    return f"{x * 100:.0f}%"


def _mechanism_paragraph(row: pd.Series) -> str:
    bits = []
    mech = row.get("mechanistic_support")
    if pd.notna(mech) and mech >= 0.3:
        bits.append(f"Drug-target convergence on the disease scores "
                    f"<b>{mech:.2f}</b> (noisy-OR over OT target-disease "
                    f"associations).")
    nov_status = _s(row.get("novelty_status"))
    if nov_status and nov_status != "novel":
        bits.append(f"Novelty is graded <i>{nov_status}</i> -- close to an "
                    f"existing indication, with the same molecular argument.")
    elif row.get("novelty", 0) >= 1.0:
        bits.append("This pair is fully novel relative to all known "
                    "drug-indication assignments in ChEMBL.")
    direction = _s(row.get("direction_status"))
    if direction == "aligned":
        bits.append("Therapeutic direction is <b>aligned</b> with OT genetic "
                    "evidence (the drug pushes the target the disease-protective "
                    "way).")
    elif direction == "opposed":
        bits.append("Therapeutic direction looks <b>opposed</b> -- a red flag "
                    "to discuss.")
    tissue = _s(row.get("tissue_status"))
    if tissue == "expressed":
        bits.append(f"Lead target is well expressed in the disease-relevant "
                    f"tissue ({row.get('tissue_evidence', '?')}).")
    phylo = row.get("phylo_score", 0) or 0
    if phylo and phylo > 0.5:
        bits.append(f"Orthologous-gene model-organism evidence is strong "
                    f"(phylo_score {phylo:.2f}, n={_fmt_int(row.get('phylo_n_models'))}).")
    pw = row.get("n_pathway_overlap", 0) or 0
    if pw:
        bits.append(f"Reactome co-membership reinforces the indirect "
                    f"mechanism with <b>{int(pw)}</b> shared specific pathways.")
    severity = _s(row.get("severity_concern"))
    if severity:
        bits.append("Note: receptor-LoF severity flag is present; mechanism "
                    "may require partial-residual-function caveat.")
    if not bits:
        bits.append("Mechanism evidence is thin -- worth probing.")
    return " ".join(bits)


def _ip_paragraph(row: pd.Series) -> str:
    loe = row.get("loe_year")
    if pd.isna(loe):
        return ("The substance is a biologic or otherwise outside the FDA "
                "Orange Book; check the Purple Book and direct registry "
                "for IP runway.")
    parts = []
    has_gen = row.get("has_generic")
    if has_gen is True or has_gen == "True":
        parts.append(f"Loss-of-exclusivity year <b>{int(loe)}</b>, with "
                     f"{_fmt_int(row.get('n_anda'))} generic ANDAs already "
                     f"on the market.")
    else:
        parts.append(f"Loss-of-exclusivity year <b>{int(loe)}</b> with "
                     f"<b>no generic competition</b> -- still meaningful IP "
                     f"runway for a co-development conversation.")
    inv = row.get("investigation_prior")
    if pd.notna(inv):
        if inv < 0.2:
            parts.append("Literature footprint is thin -- a first-mover "
                         "publication is plausible.")
        elif inv > 0.7:
            parts.append("Literature footprint is dense -- competition is "
                         "likely; positioning matters.")
    return " ".join(parts)


def _market_paragraph(row: pd.Series) -> str:
    pts = row.get("us_patients")
    if pd.isna(pts):
        return ("US prevalence figures are not in the curated map; for an "
                "orphan condition this is expected. The KOL can size the "
                "addressable population.")
    is_orph = row.get("is_orphan")
    src = _s(row.get("market_source"))
    if is_orph is True or is_orph == "True":
        return (f"US population: approximately <b>{_fmt_int(pts)}</b> "
                f"patients -- below the 200,000 FDA Orphan Drug Act "
                f"threshold. Orphan exclusivity + premium pricing changes "
                f"the deal math materially. (source: {src})")
    return (f"US population: approximately <b>{_fmt_int(pts)}</b> "
            f"patients -- a mass-market opportunity. (source: {src})")


def _kol_paragraph(row: pd.Series, region: str) -> tuple[str, str, str]:
    """(salutation, body, signature line) for one region. Returns ('','','')
    if no KOL was identified for that region."""
    name = _s(row.get(f"{region}_kol_name"))
    if not name:
        return "", "", ""
    inst = _s(row.get(f"{region}_kol_institution"))
    email = _s(row.get(f"{region}_kol_email"))
    h = row.get(f"{region}_kol_h_index")
    n_pubs = row.get(f"{region}_kol_n_pubs")

    salut = f"For the attention of: {name}"
    if inst:
        salut += f", {inst}"
    creds = []
    if pd.notna(h):
        creds.append(f"h-index {int(h)}")
    if pd.notna(n_pubs):
        creds.append(f"{int(n_pubs)} relevant PubMed publication{'s' if int(n_pubs) > 1 else ''}")
    body = (f"You were identified from a PubMed pass over recent papers "
            f"co-mentioning the target and the indication"
            + (f" ({'; '.join(creds)})" if creds else "") + ".")
    contact = email if email else ""
    return salut, body, contact


# ----------------------------------------------------------------------
# PDF rendering
# ----------------------------------------------------------------------
def _safe_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s.strip())[:80] or "memo"


def _render_memo(row: pd.Series, region: str, profile: dict, out_dir: Path,
                 *, no_kol: bool = False) -> Path | None:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.lib.enums import TA_LEFT, TA_CENTER
    from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                    HRFlowable)
    from reportlab.lib import colors

    if no_kol:
        salut, kol_body, kol_contact = "For clinical-collaborator outreach.", "", ""
    else:
        salut, kol_body, kol_contact = _kol_paragraph(row, region)
        if not salut:
            return None  # Skip regions without a KOL

    substance = _s(row.get("substance_name")) or _s(row.get("drug_name")) or "?"
    disease = _s(row.get("disease_name")) or "?"
    safe = f"{_safe_filename(substance)}__{_safe_filename(disease)}__{region}.pdf"
    out_path = out_dir / safe

    styles = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=styles["Normal"], fontSize=9,
                          leading=12, spaceAfter=4, alignment=TA_LEFT)
    section = ParagraphStyle("section", parent=styles["Heading4"], fontSize=10,
                             leading=13, textColor=colors.HexColor("#0B4FA8"),
                             spaceBefore=8, spaceAfter=2)
    title = ParagraphStyle("title", parent=styles["Title"], fontSize=15,
                           leading=18, alignment=TA_LEFT, spaceAfter=4)
    meta = ParagraphStyle("meta", parent=styles["Normal"], fontSize=8,
                          leading=10, textColor=colors.HexColor("#666666"))
    foot = ParagraphStyle("foot", parent=styles["Normal"], fontSize=7,
                          leading=9, textColor=colors.HexColor("#888888"),
                          alignment=TA_CENTER)

    doc = SimpleDocTemplate(str(out_path), pagesize=A4,
                            leftMargin=2 * cm, rightMargin=2 * cm,
                            topMargin=1.5 * cm, bottomMargin=1.5 * cm,
                            title=f"{substance} -> {disease}")
    story = []
    story.append(Paragraph("REPURPOSING COLLABORATION OPPORTUNITY", meta))
    story.append(Paragraph(f"<b>{substance}</b> &nbsp;&nbsp;->&nbsp;&nbsp; "
                           f"<b>{disease}</b>", title))
    story.append(Paragraph(salut, body))
    story.append(HRFlowable(width="100%", thickness=0.6,
                            color=colors.HexColor("#cccccc"), spaceAfter=4))

    story.append(Paragraph("Why this hypothesis is real", section))
    story.append(Paragraph(_mechanism_paragraph(row), body))

    story.append(Paragraph("Why now (IP runway + literature footprint)", section))
    story.append(Paragraph(_ip_paragraph(row), body))

    story.append(Paragraph("Clinical opportunity", section))
    story.append(Paragraph(_market_paragraph(row), body))

    if kol_body:
        story.append(Paragraph("Why you", section))
        story.append(Paragraph(kol_body, body))
        if kol_contact:
            story.append(Paragraph(f"On record: <font face='Courier'>{kol_contact}</font>", body))

    story.append(Paragraph("What we propose", section))
    g = profile["group"]
    pitch_lines = "<br/>".join(f"&bull; {line}" for line in g.get("pitch", []))
    story.append(Paragraph(
        f"<b>{g['name']}</b> ({g['affiliation']}) brings bioinformatics "
        f"support to the collaboration:<br/>{pitch_lines}", body))
    if profile.get("collaboration_model"):
        story.append(Paragraph("Proposed collaboration model:", body))
        model_lines = "<br/>".join(f"&bull; {line}"
                                    for line in profile["collaboration_model"])
        story.append(Paragraph(model_lines, body))

    story.append(Paragraph("Next step", section))
    story.append(Paragraph(profile["next_step"]["ask"].strip(), body))

    story.append(Spacer(1, 0.4 * cm))
    sig = (f"<b>{g.get('pi_name', '')}</b>"
           + (f", {g['pi_title']}" if g.get("pi_title") else "")
           + (f" ({g['pi_credentials']})" if g.get("pi_credentials") else "")
           + f"<br/>{g.get('affiliation', '')}"
           + (f"<br/><font face='Courier'>{g['pi_email']}</font>"
              if g.get("pi_email") else ""))
    story.append(Paragraph(sig, body))

    story.append(Spacer(1, 0.4 * cm))
    story.append(HRFlowable(width="100%", thickness=0.4,
                            color=colors.HexColor("#dddddd"), spaceAfter=4))
    story.append(Paragraph(profile["footer_disclaimer"].strip()
                           .replace("\n", " "), foot))
    story.append(Paragraph(f"Generated {date.today().isoformat()}.", foot))

    doc.build(story)
    return out_path


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DEFAULT_INPUT))
    ap.add_argument("--output", default=str(DEFAULT_OUTPUT))
    ap.add_argument("--profile", default=str(DEFAULT_PROFILE))
    ap.add_argument("--top", type=int, default=5,
                    help="How many hypotheses (with at least one KOL) to render")
    ap.add_argument("--region", choices=["us", "eu", "both"], default="both")
    ap.add_argument("--orphan-only", action="store_true")
    ap.add_argument("--filter-min-opportunity", type=float, default=0.0)
    ap.add_argument("--pi-name", default=None)
    ap.add_argument("--pi-email", default=None)
    ap.add_argument("--group-name", default=None)
    ap.add_argument("--no-kol", action="store_true",
                    help="Skip KOL identification; emit one memo per top hypothesis "
                         "with a generic salutation suitable for blinded outreach.")
    args = ap.parse_args()

    profile = _load_profile(Path(args.profile), {
        "pi_name": args.pi_name, "pi_email": args.pi_email,
        "group_name": args.group_name,
    })
    df = pd.read_csv(args.input, low_memory=False)
    if args.filter_min_opportunity:
        df = df[df["opportunity"] >= args.filter_min_opportunity]
    if args.orphan_only:
        df = df[df["is_orphan"] == True]
    if not args.no_kol:
        # We need at least one KOL on the selected side(s).
        if args.region == "us":
            df = df[df["us_kol_name"].fillna("") != ""]
        elif args.region == "eu":
            df = df[df["eu_kol_name"].fillna("") != ""]
        else:
            df = df[(df["us_kol_name"].fillna("") != "") |
                    (df["eu_kol_name"].fillna("") != "")]
    df = df.sort_values("opportunity", ascending=False).head(args.top)

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    regions = ["us", "eu"] if args.region == "both" else [args.region]
    if args.no_kol:
        regions = regions[:1]  # one memo per hypothesis, region-agnostic
    written: list[Path] = []
    for _, row in df.iterrows():
        for region in regions:
            p = _render_memo(row, region, profile, out_dir, no_kol=args.no_kol)
            if p is not None:
                written.append(p)

    print(f"wrote {len(written)} memos to {out_dir}")
    for p in written:
        try:
            print(f"  {p.resolve().relative_to(ROOT)}")
        except ValueError:
            print(f"  {p}")


if __name__ == "__main__":
    main()
