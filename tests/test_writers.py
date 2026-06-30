"""Tests for the cultivar coefficient writer (no DSSAT execution needed)."""
import shutil
from pathlib import Path

from dssatcalibrator.writers import edit_cultivar, edit_species, read_cultivar_values, cultivar_field_map

GENO = Path("C:/DSSAT48/Genotype/HMGRO048.CUL")


def _cul(tmp_path):
    import pytest
    if not GENO.exists():
        pytest.skip("HMGRO048.CUL not present")
    dst = tmp_path / "HMGRO048.CUL"
    shutil.copy(GENO, dst)
    return dst


def test_field_map_has_all_coefficients(tmp_path):
    cul = _cul(tmp_path)
    fmap = cultivar_field_map(cul)
    assert {"CSDL", "PPSEN", "EM-FL", "FL-SD", "SD-PM", "LFMAX", "SLAVR", "SIZLF"} <= set(fmap)
    # every field is the same fixed width DSSAT uses
    widths = {hi - lo for lo, hi in fmap.values()}
    assert widths == {6}


def test_edit_cultivar_changes_only_targeted(tmp_path):
    cul = _cul(tmp_path)
    before = read_cultivar_values(cul, "IB0008")
    edit_cultivar(cul, "IB0008", {"CSDL": 14.5, "EM-FL": 18.0})
    after = read_cultivar_values(cul, "IB0008")
    assert after["CSDL"] == 14.5
    assert after["EM-FL"] == 18.0
    # untouched coefficients preserved exactly
    for k in ("PPSEN", "LFMAX", "SLAVR", "SIZLF", "WTPSD", "THRSH"):
        assert after[k] == before[k], k


def test_edit_cultivar_rejects_unknown(tmp_path):
    import pytest
    cul = _cul(tmp_path)
    with pytest.raises(KeyError):
        edit_cultivar(cul, "IB0008", {"NOPE": 1.0})


def test_edit_species_preserves_leading_dot_decimal(tmp_path):
    spe = tmp_path / "sample.SPE"
    spe.write_text(
        " .0046 .0004 .3000  4.90 1.030             SLWREF,SLWSLO,NSLOPE,LNREF,PGREF\n"
        "! .0036 .0004 .3000  4.90 1.030             SLWREF,SLWSLO,NSLOPE,LNREF,PGREF\n",
        encoding="utf-8",
    )

    edit_species(spe, {"SLWREF,SLWSLO,NSLOPE,LNREF,PGREF": 0.0048})

    lines = spe.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith(" .0048 .0004")
    assert lines[1].startswith("! .0036 .0004")


def test_edit_species_can_target_numeric_token_index(tmp_path):
    spe = tmp_path / "sample.SPE"
    spe.write_text(
        "  40.00 61.00  0.96  0.10                   PARMAX,PHTMAX,KCAN,KC_SLOPE\n",
        encoding="utf-8",
    )

    edit_species(spe, {"PARMAX,PHTMAX,KCAN,KC_SLOPE": {"value": 0.88, "index": 2}})

    line = spe.read_text(encoding="utf-8").splitlines()[0]
    assert line.startswith("  40.00 61.00  0.88  0.10")
