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
