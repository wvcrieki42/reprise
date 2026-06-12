"""Unit tests for FDA Orange Book ingestion."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from repurpose.sources import adapters  # noqa: E402
from repurpose.sources.loaders import load_orange_book  # noqa: E402


def test_normalize_ingredient_strips_salt():
    assert adapters._normalize_ingredient("PIOGLITAZONE HYDROCHLORIDE") == "PIOGLITAZONE"
    assert adapters._normalize_ingredient("METFORMIN HYDROCHLORIDE") == "METFORMIN"
    # Stacked suffixes ("X SULFATE MONOHYDRATE") get peeled one at a time
    assert adapters._normalize_ingredient("DRUG SULFATE MONOHYDRATE") == "DRUG"
    # No suffix -> unchanged (uppercase + strip)
    assert adapters._normalize_ingredient(" lubiprostone ") == "LUBIPROSTONE"
    # Empty / None
    assert adapters._normalize_ingredient("") == ""
    assert adapters._normalize_ingredient(None) == ""


def _make_ob_dir(tmp_path: Path) -> Path:
    """Synthesize tiny Orange Book TSV files mirroring the FDA schema."""
    d = tmp_path / "ob"
    d.mkdir()
    # Two NDAs (brand) + one ANDA (generic) for PIOGLITAZONE HYDROCHLORIDE,
    # one NDA for VOXELOTOR (no generics), one combination product.
    (d / "products.txt").write_text(
        "Ingredient~DF;Route~Trade_Name~Applicant~Strength~Appl_Type~Appl_No~Product_No~TE_Code~Approval_Date~RLD~RS~Type~Applicant_Full_Name\n"
        "PIOGLITAZONE HYDROCHLORIDE~TABLET;ORAL~ACTOS~TAKEDA~15MG~N~021073~001~AB~Jul 1999~Yes~Yes~RX~TAKEDA\n"
        "PIOGLITAZONE HYDROCHLORIDE~TABLET;ORAL~PIOGLITAZONE~MYLAN~15MG~A~090001~001~AB~Aug 2012~No~No~RX~MYLAN\n"
        "VOXELOTOR~TABLET;ORAL~OXBRYTA~PFIZER~500MG~N~213137~001~~Nov 2019~Yes~Yes~RX~PFIZER\n"
        "METFORMIN HYDROCHLORIDE; SITAGLIPTIN PHOSPHATE~TABLET;ORAL~JANUMET~MERCK~50MG/500MG~N~022044~001~AB~Mar 2007~Yes~Yes~RX~MERCK\n"
    )
    (d / "patent.txt").write_text(
        "Appl_Type~Appl_No~Product_No~Patent_No~Patent_Expire_Date_Text~Drug_Substance_Flag~Drug_Product_Flag~Patent_Use_Code~Delist_Flag~Submission_Date\n"
        "N~021073~001~5965584~Aug 24, 2026~~~~~\n"
        "N~021073~001~6291490~Sep 12, 2029~~~~~\n"
        "N~213137~001~9518056~Mar 15, 2037~~~~~\n"
        "N~022044~001~7625884~Dec 12, 2026~~~~~\n"
    )
    (d / "exclusivity.txt").write_text(
        "Appl_Type~Appl_No~Product_No~Exclusivity_Code~Exclusivity_Date\n"
        "N~213137~001~ODE~Nov 15, 2026\n"
    )
    return d


def test_adapter_per_ingredient_aggregation(tmp_path):
    ob_dir = _make_ob_dir(tmp_path)
    df = adapters.fda_orange_book(str(ob_dir))

    pio = df[df.ingredient == "PIOGLITAZONE"].iloc[0]
    # Two patents: 2026 and 2029 -> latest = 2029
    assert int(pio["latest_patent_year"]) == 2029
    assert pd.isna(pio["latest_exclusivity_year"])
    assert int(pio["loe_year"]) == 2029
    assert pio["has_generic"]                # 1 ANDA exists
    assert int(pio["n_nda"]) == 1
    assert int(pio["n_anda"]) == 1

    vox = df[df.ingredient == "VOXELOTOR"].iloc[0]
    assert int(vox["latest_patent_year"]) == 2037
    assert int(vox["latest_exclusivity_year"]) == 2026
    assert int(vox["loe_year"]) == 2037     # max(patent, excl)
    assert not vox["has_generic"]            # NDA only, no ANDA
    assert int(vox["n_anda"]) == 0


def test_adapter_explodes_combinations(tmp_path):
    """A combination product 'METFORMIN; SITAGLIPTIN' should contribute to BOTH ingredients."""
    ob_dir = _make_ob_dir(tmp_path)
    df = adapters.fda_orange_book(str(ob_dir))
    # Both ingredient names appear in the output (normalised, no salt suffix)
    assert "METFORMIN" in df.ingredient.tolist()
    assert "SITAGLIPTIN" in df.ingredient.tolist()
    # The same patent (Dec 2026) is attributed to both
    metformin = df[df.ingredient == "METFORMIN"].iloc[0]
    sita = df[df.ingredient == "SITAGLIPTIN"].iloc[0]
    assert int(metformin["latest_patent_year"]) == 2026
    assert int(sita["latest_patent_year"]) == 2026


def test_loader_missing_file_returns_empty(tmp_path):
    out = load_orange_book(tmp_path / "missing.csv")
    assert out.empty
    assert "loe_year" in out.columns


def test_loader_round_trip(tmp_path):
    csv = tmp_path / "ob.csv"
    csv.write_text(
        "ingredient,latest_patent_year,latest_exclusivity_year,loe_year,has_generic,n_nda,n_anda\n"
        "PIOGLITAZONE,2029,,2029,True,1,1\n"
        "VOXELOTOR,2037,2026,2037,False,3,0\n"
    )
    out = load_orange_book(csv)
    pio = out[out.ingredient == "PIOGLITAZONE"].iloc[0]
    assert int(pio["latest_patent_year"]) == 2029
    assert pio["has_generic"] is True or pio["has_generic"] == True
    vox = out[out.ingredient == "VOXELOTOR"].iloc[0]
    assert vox["has_generic"] is False or vox["has_generic"] == False


if __name__ == "__main__":
    import tempfile
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            sig = fn.__code__.co_varnames[:fn.__code__.co_argcount]
            with tempfile.TemporaryDirectory() as d:
                if "tmp_path" in sig:
                    fn(Path(d))
                else:
                    fn()
                print(f"OK {name}")
