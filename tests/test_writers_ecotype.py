"""Tests for the ecotype coefficient writer (no DSSAT execution needed)."""
import shutil
from pathlib import Path

from dssatcalibrator.writers import edit_ecotype, read_ecotype_values, ecotype_field_map

ECO = Path("C:/DSSAT48/Genotype/HMGRO048.ECO")


def _eco(tmp_path):
    import pytest
    if not ECO.exists():
        pytest.skip("HMGRO048.ECO not present")
    dst = tmp_path / "HMGRO048.ECO"
    shutil.copy(ECO, dst)
    return dst


def test_field_map_has_all_coefficients(tmp_path):
    eco = _eco(tmp_path)
    fmap = ecotype_field_map(eco)
    assert {"PL-EM", "RHGHT", "EM-V1", "FL-VS", "JU-R0", "LNGSH", "TRIFL"} <= set(fmap)
    # every field is the same fixed width DSSAT uses
    widths = {hi - lo for lo, hi in fmap.values()}
    assert widths == {3, 6}


def test_edit_ecotype_changes_only_targeted(tmp_path):
    eco = _eco(tmp_path)
    before = read_ecotype_values(eco, "HM0003")
    edit_ecotype(eco, "HM0003", {"PL-EM": 5.2, "RHGHT": 2.1})
    after = read_ecotype_values(eco, "HM0003")
    assert after["PL-EM"] == 5.2
    assert after["RHGHT"] == 2.1
    # untouched coefficients preserved exactly
    for k in ("EM-V1", "FL-VS", "JU-R0", "LNGSH", "TRIFL"):
        assert after[k] == before[k], k


def test_edit_ecotype_rejects_unknown(tmp_path):
    import pytest
    eco = _eco(tmp_path)
    with pytest.raises(KeyError):
        edit_ecotype(eco, "HM0003", {"NOPE": 1.0})


def test_edit_ecotype_keeps_precision_for_small_changed_values(tmp_path):
    eco = tmp_path / "HMGRO048.ECO"
    eco.write_text(
        "*HEMP ECOTYPE COEFFICIENTS: CRGRO048 MODEL\n"
        "@ECO#  ECONAME.......... MG TM THVAR PL-EM EM-V1 V1-JU JU-R0  PM06  PM09 LNGSH R7-R8 FL-VS TRIFL RWDTH RHGHT R1PPO OPTBI SLOBI\n"
        "HM0003 HYBRID HEMP China 07 01   0.0   4.0   3.0   0.0   1.0   0.0  0.80   5.0  20.0 50.00  0.25   1.0   1.5  .504  20.0  .035\n",
        encoding="utf-8",
    )

    edit_ecotype(eco, "HM0003", {"THVAR": 0.0393, "SLOBI": 0.03442})
    after = read_ecotype_values(eco, "HM0003")

    assert after["THVAR"] == 0.0393
    assert after["SLOBI"] == 0.0344

    row = next(ln for ln in eco.read_text(encoding="utf-8").splitlines() if ln.startswith("HM0003"))
    fmap = ecotype_field_map(eco)
    assert row[fmap["SLOBI"][0]:fmap["SLOBI"][1]].startswith(" ")
