"""End-to-end smoke test on the bundled sample data.

Run:  python -m pytest -q   (or)   python tests/test_smoke.py
Checks core logic: universe filtering, novelty subtraction, network expansion,
the directionality model (aligned vs opposed), and pandas/duckdb engine parity.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.config import load_config          # noqa: E402
from repurpose.pipeline import run                 # noqa: E402


def _run(overrides=None):
    cfg = load_config(ROOT / "config.yaml")
    if overrides:
        for path, val in overrides.items():
            node = cfg.raw
            *parents, leaf = path.split(".")
            for p in parents:
                node = node.setdefault(p, {})
            node[leaf] = val
    return run(cfg, verbose=False)


def test_core_logic():
    df = _run()
    assert len(df) > 0
    assert {"drug_name", "disease_name", "opportunity", "novelty", "direction_factor",
            "tissue_factor", "lead_target", "lead_target_name",
            "lead_target_genecards", "data_version"} <= set(df.columns)

    # lead-target name + GeneCards link populated
    a_ipf = df[(df.drug_name == "Auranofin") & (df.disease_name == "Idiopathic pulmonary fibrosis")].iloc[0]
    assert a_ipf["lead_target"] == "TXNRD1"
    assert "reductase" in a_ipf["lead_target_name"].lower()
    assert a_ipf["lead_target_genecards"].endswith("gene=TXNRD1")
    # provenance recorded on every row
    assert df["data_version"].str.contains("ChEMBL").all()

    # 1. non-approved drug excluded
    assert "Experimental-NotApproved" not in set(df["drug_name"])

    # 2. known indications subtracted
    pairs = set(zip(df["drug_name"], df["disease_name"]))
    assert ("Auranofin", "Rheumatoid arthritis") not in pairs
    assert ("Sildenafil", "Pulmonary arterial hypertension") not in pairs

    # 3. network expansion surfaced neighbour-only diseases for Auranofin
    aura = df[df["drug_name"] == "Auranofin"]
    aura_dis = set(aura["disease_name"])
    assert {"Amyotrophic lateral sclerosis", "Ulcerative colitis",
            "Idiopathic pulmonary fibrosis"} <= aura_dis

    # 4. documented lead still tops Auranofin
    top = aura.sort_values("opportunity", ascending=False).iloc[0]["disease_name"]
    assert top == "Idiopathic pulmonary fibrosis", f"expected IPF, got {top}"
    print(f"OK core: {len(df)} hypotheses; Auranofin lead = {top}")


def test_directionality():
    df = _run()
    # aligned: Auranofin->IPF (inhibitor of targets where inhibition is therapeutic)
    a_ipf = df[(df.drug_name == "Auranofin") & (df.disease_name == "Idiopathic pulmonary fibrosis")].iloc[0]
    assert a_ipf["direction_status"] == "aligned"
    assert a_ipf["direction_factor"] == 1.0

    # opposed: Metformin ACTIVATES AMPK where the (illustrative) therapeutic direction is to lower it
    m_ad = df[(df.drug_name == "Metformin") & (df.disease_name == "Alzheimer disease")]
    assert len(m_ad) == 1
    assert m_ad.iloc[0]["direction_status"] == "opposed"
    assert m_ad.iloc[0]["direction_factor"] < 0.2

    # disabling direction removes the penalty -> Metformin->AD opportunity rises
    on = df[(df.drug_name == "Metformin") & (df.disease_name == "Alzheimer disease")].iloc[0]["opportunity"]
    df_off = _run({"direction.enabled": False})
    off = df_off[(df_off.drug_name == "Metformin") & (df_off.disease_name == "Alzheimer disease")].iloc[0]["opportunity"]
    assert off > on, "disabling directionality should raise the opposed hypothesis's score"
    print(f"OK direction: opposed penalised ({on:.4f}) vs neutral ({off:.4f})")


def test_tissue_filter():
    df = _run()
    # expressed: Auranofin->IPF targets are present in lung
    a_ipf = df[(df.drug_name == "Auranofin") & (df.disease_name == "Idiopathic pulmonary fibrosis")].iloc[0]
    assert a_ipf["tissue_status"] == "expressed"
    assert a_ipf["tissue_factor"] == 1.0

    # penalised: Simvastatin->Alzheimer relies on HMGCR, which is barely expressed in brain
    s_ad = df[(df.drug_name == "Simvastatin") & (df.disease_name == "Alzheimer disease")]
    assert len(s_ad) == 1
    assert s_ad.iloc[0]["tissue_status"] in {"low", "absent"}
    assert s_ad.iloc[0]["tissue_factor"] < 1.0

    # turning the tissue gate off raises that penalised hypothesis
    on = s_ad.iloc[0]["opportunity"]
    off = _run({"tissue.enabled": False})
    off_v = off[(off.drug_name == "Simvastatin") & (off.disease_name == "Alzheimer disease")].iloc[0]["opportunity"]
    assert off_v > on
    print(f"OK tissue: HMGCR/brain penalised ({on:.4f}) vs off ({off_v:.4f})")


def test_engine_parity():
    pd_df = _run({"engine": "pandas"})
    try:
        dk_df = _run({"engine": "duckdb"})
    except ImportError:
        print("SKIP parity: duckdb not installed")
        return
    pk = set(zip(pd_df.drug_name, pd_df.disease_name))
    dk = set(zip(dk_df.drug_name, dk_df.disease_name))
    assert pk == dk, f"engines disagree on hypothesis set: {pk ^ dk}"
    # opportunity values agree within tolerance
    merged = pd_df.merge(dk_df, on=["drug_id", "efo_id"], suffixes=("_pd", "_dk"))
    assert len(merged) == len(pd_df)
    diff = (merged["opportunity_pd"] - merged["opportunity_dk"]).abs().max()
    assert diff < 1e-4, f"opportunity mismatch up to {diff}"
    print(f"OK parity: {len(pd_df)} hypotheses identical, max |Δopportunity| = {diff:.2e}")


if __name__ == "__main__":
    test_core_logic()
    test_directionality()
    test_tissue_filter()
    test_engine_parity()
    print("all smoke tests passed")
