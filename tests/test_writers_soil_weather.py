"""Tests for the soil (.SOL) and weather (.WTH) writers + FileX field parsing."""
import shutil
from pathlib import Path

import pytest

from dssatcalibrator.writers import (
    edit_soil,
    edit_weather,
    extract_soil_profile,
    parse_fields,
    parse_header_boundaries,
)

REPO = Path(__file__).resolve().parents[1]
SMOKE_HMX = REPO / "_smoke" / "YUKU2101.HMX"
DSSAT = Path("C:/DSSAT48")
SOIL_SOL = DSSAT / "Soil" / "SOIL.SOL"
WTH = DSSAT / "Weather" / "CNKU2101.WTH"


def test_parse_fields():
    if not SMOKE_HMX.exists():
        pytest.skip("YUKU2101.HMX not present in _smoke")
    f = parse_fields(SMOKE_HMX)
    assert f["id_soil"] == "YUKU2101"
    assert f["wsta"] == "CNKU2101"


def test_parse_fields_coordinates(tmp_path):
    filex = tmp_path / "TEST.HMX"
    filex.write_text(
        "*FIELDS\n"
        "@L ID_FIELD WSTA....  FLSA  FLOB  FLDT  FLDD  FLDS  FLST SLTX  SLDP  ID_SOIL    FLNAME\n"
        " 1 FIELD001 WSTA0001   -99     0 IB000     0     0 00000 -99    180  SOIL0001   test\n"
        "@L ...........XCRD ...........YCRD .....ELEV .............AREA .SLEN .FLWR .SLAS FLHST FHDUR\n"
        " 1        -91.250        14.500      120               -99   -99   -99   -99   -99   -99\n"
        "*TREATMENTS\n",
        encoding="utf-8",
    )
    f = parse_fields(filex)
    assert f["id_field"] == "FIELD001"
    assert f["wsta"] == "WSTA0001"
    assert f["id_soil"] == "SOIL0001"
    assert f["lat"] == pytest.approx(14.5)
    assert f["lon"] == pytest.approx(-91.25)
    assert f["elev"] == pytest.approx(120)


def _layer_values(sol_text: str, col: str):
    lines = sol_text.splitlines()
    hdr = next(i for i, ln in enumerate(lines)
               if ln.lstrip().startswith("@") and "SDUL" in ln and "SLB" in ln)
    fmap = parse_header_boundaries(lines[hdr])
    lo, hi = fmap[col]
    out = []
    for ln in lines[hdr + 1:]:
        if ln.lstrip().startswith(("@", "*")):
            break
        if ln.strip() and ln.lstrip()[0].isdigit():
            out.append(float(ln[lo:hi]))
    return out


def _profile_value(sol_text: str, col: str):
    lines = sol_text.splitlines()
    hdr = next(i for i, ln in enumerate(lines)
               if ln.lstrip().startswith("@") and "SALB" in ln and "SDUL" not in ln)
    fmap = parse_header_boundaries(lines[hdr])
    lo, hi = fmap[col]
    return float(lines[hdr + 1][lo:hi])


def test_extract_and_edit_soil(tmp_path):
    if not SOIL_SOL.exists():
        pytest.skip("SOIL.SOL not present")
    try:
        block = extract_soil_profile(SOIL_SOL, "YUKU2101")
    except ValueError as exc:
        pytest.skip(str(exc))
    assert block.startswith("*YUKU2101")
    # only one profile in the extracted block
    assert block.count("\n*") == 0

    dst = tmp_path / "SOIL.SOL"
    dst.write_text(block)

    sdul_before = _layer_values(dst.read_text(), "SDUL")
    edit_soil(dst, "YUKU2101", layer_mults={"SDUL": 1.1}, profile_sets={"SLPF": 0.8})
    after = dst.read_text()
    sdul_after = _layer_values(after, "SDUL")

    assert len(sdul_after) == len(sdul_before) == 3
    for b, a in zip(sdul_before, sdul_after):
        assert a == pytest.approx(min(b * 1.1, 1.0), abs=1e-3)
    assert _profile_value(after, "SLPF") == pytest.approx(0.8, abs=1e-3)


def test_edit_weather(tmp_path):
    if not WTH.exists():
        pytest.skip("CNKU2101.WTH not present")
    dst = tmp_path / "CNKU2101.WTH"
    shutil.copy(WTH, dst)

    def first_data_row(text):
        # weather data rows start with a 7-digit YYYYDDD date
        for ln in text.splitlines():
            if ln[:7].strip().isdigit() and len(ln[:7].strip()) == 7:
                return ln
        return None

    fmap = parse_header_boundaries(
        next(ln for ln in dst.read_text().splitlines()
             if ln.lstrip().startswith("@") and "SRAD" in ln))
    slo, shi = fmap["SRAD"]
    tlo, thi = fmap["TMAX"]
    before = first_data_row(dst.read_text())
    srad_before = float(before[slo:shi])
    tmax_before = float(before[tlo:thi])

    edit_weather(dst, {"SRAD": ("mult", 1.1), "TMAX": ("off", 1.0)})
    after = first_data_row(dst.read_text())
    assert after is not None
    # columns remain whitespace-separated (free-format safe)
    assert len(after.split()) >= 5
    assert float(after[slo:shi]) == pytest.approx(srad_before * 1.1, abs=0.06)
    assert float(after[tlo:thi]) == pytest.approx(tmax_before + 1.0, abs=0.06)
