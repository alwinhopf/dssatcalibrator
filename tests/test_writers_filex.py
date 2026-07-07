"""Tests for the FileX management and initial conditions writers."""
import shutil
import re
from pathlib import Path

import pytest

from dssatcalibrator.writers import edit_filex, parse_header_boundaries

REPO = Path(__file__).resolve().parents[1]
SMOKE_HMX = REPO / "_smoke" / "YUKU2101.HMX"


def _row(header, values):
    fmap = parse_header_boundaries(header)
    width = max(hi for _, hi in fmap.values())
    chars = [" "] * width
    for field, value in values.items():
        lo, hi = fmap[field]
        text = str(value).strip().rjust(hi - lo)
        if len(text) > hi - lo:
            raise ValueError(f"{field} value {value!r} too wide for test row")
        chars[lo:hi] = list(text)
    return "".join(chars)


def _cell(lines, section, header_prefix, field, row=1):
    sec = next(i for i, ln in enumerate(lines) if ln.startswith(section))
    hdr = next(i for i in range(sec + 1, len(lines)) if lines[i].lstrip().startswith(header_prefix))
    fmap = parse_header_boundaries(lines[hdr])
    seen = 0
    for i in range(hdr + 1, len(lines)):
        ln = lines[i]
        if ln.startswith("*") or ln.lstrip().startswith("@"):
            break
        if not ln.strip() or ln.lstrip().startswith("!"):
            continue
        seen += 1
        if seen == row:
            lo, hi = fmap[field]
            return ln[lo:hi].strip()
    raise AssertionError(f"row {row} not found in {section}")


def _synthetic_filex(path):
    p_hdr = "@P PDATE EDATE  PPOP  PPOE  PLME  PLDS  PLRS  PLRD  PLDP"
    c_hdr = "@C  PCR ICDAT  ICBL  SH2O  SNH4  SNO3"
    i_hdr = "@I IDATE IROP  IRVAL"
    f_hdr = "@F FDATE FMCD FACD FDEP FAMN FAMP FAMK"
    path.write_text(
        "\n".join([
            "*SYNTHETIC FILEX",
            "*TREATMENTS",
            "@N R O C TNAME....................",
            " 1 1 0 0 CONTROL",
            "*PLANTING DETAILS",
            p_hdr,
            _row(p_hdr, {"P": 1, "PDATE": 21150, "EDATE": -99, "PPOP": 30.0, "PPOE": 30.0,
                         "PLME": "S", "PLDS": "R", "PLRS": 50.0, "PLRD": -99, "PLDP": 3.0}),
            _row(p_hdr, {"P": 2, "PDATE": 21151, "EDATE": -99, "PPOP": 31.0, "PPOE": 31.0,
                         "PLME": "S", "PLDS": "R", "PLRS": 55.0, "PLRD": -99, "PLDP": 3.5}),
            "*INITIAL CONDITIONS",
            c_hdr,
            _row(c_hdr, {"C": 1, "PCR": 100, "ICDAT": 21140, "ICBL": 5,
                         "SH2O": ".200", "SNH4": ".010", "SNO3": 1.20}),
            _row(c_hdr, {"C": 1, "PCR": 100, "ICDAT": 21140, "ICBL": 30,
                         "SH2O": ".250", "SNH4": ".020", "SNO3": 1.40}),
            "*IRRIGATION AND WATER MANAGEMENT",
            i_hdr,
            _row(i_hdr, {"I": 1, "IDATE": 21155, "IROP": "IR001", "IRVAL": 10.0}),
            _row(i_hdr, {"I": 1, "IDATE": 21170, "IROP": "IR001", "IRVAL": 15.0}),
            "*FERTILIZERS (INORGANIC)",
            f_hdr,
            _row(f_hdr, {"F": 1, "FDATE": 21156, "FMCD": "FE001", "FACD": "AP001",
                         "FDEP": 5.0, "FAMN": 30.0, "FAMP": 0.0, "FAMK": 0.0}),
        ]) + "\n",
        encoding="utf-8",
    )


def test_parse_header_boundaries():
    header = "@P PDATE EDATE  PPOP  PPOE  PLME  PLDS  PLRS  PLRD  PLDP"
    bounds = parse_header_boundaries(header)
    
    assert "PDATE" in bounds
    assert "PPOP" in bounds
    assert "PLRS" in bounds
    assert "PLDP" in bounds
    
    # Check bounds of PPOP
    # "@P PDATE EDATE  PPOP"
    # prev token is EDATE, ends at 20
    # PPOP ends at 26
    # Exact bounds:
    assert bounds["PPOP"] == (header.index("EDATE") + 5, header.index("PPOP") + 4)
    assert parse_header_boundaries("@L ID_FIELD WSTA....  ID_SOIL")["WSTA"] == (11, 20)


def test_edit_filex_management_and_initial_water(tmp_path):
    if not SMOKE_HMX.exists():
        import pytest
        pytest.skip("YUKU2101.HMX not present in _smoke")
        
    dst = tmp_path / "YUKU2101.HMX"
    shutil.copy(SMOKE_HMX, dst)
    
    # Perturb planting details (PPOP = 45.0, PLRS = 35.0, PLDP = 4.2)
    # and initial soil water multiplier = 1.25
    edit_filex(dst, {"PPOP": 45.0, "PLRS": 35.0, "PLDP": 4.2}, {"initial_soil_water_mult": 1.25})
    
    content = dst.read_text(errors="replace")
    lines = content.splitlines()
    
    # 1. Verify planting details changes
    # Locate *PLANTING DETAILS and read the data lines
    planting_sec_idx = next(i for i, ln in enumerate(lines) if ln.startswith("*PLANTING DETAILS"))
    header_idx = next(i for i in range(planting_sec_idx + 1, len(lines)) if lines[i].lstrip().startswith("@P"))
    
    fmap = parse_header_boundaries(lines[header_idx])
    
    data_lines = []
    for i in range(header_idx + 1, len(lines)):
        ln = lines[i]
        if ln.startswith("*") or ln.lstrip().startswith("@"):
            break
        if ln.lstrip().startswith("!") or not ln.strip():
            continue
        if ln.strip() and ln.strip()[0].isdigit():
            data_lines.append(ln)
            
    assert len(data_lines) == 4  # there are 4 treatments in YUKU2101.HMX
    for ln in data_lines:
        ppop_val = float(ln[fmap["PPOP"][0]:fmap["PPOP"][1]].strip())
        plrs_val = float(ln[fmap["PLRS"][0]:fmap["PLRS"][1]].strip())
        pldp_val = float(ln[fmap["PLDP"][0]:fmap["PLDP"][1]].strip())
        
        assert ppop_val == 45.0
        assert plrs_val == 35.0
        assert pldp_val == 4.2

    # 2. Verify initial conditions (SH2O multiplied by 1.25)
    # Locate *INITIAL CONDITIONS, get @C containing SH2O
    init_sec_idx = next(i for i, ln in enumerate(lines) if ln.startswith("*INITIAL CONDITIONS"))
    c_header_idx = next(i for i in range(init_sec_idx + 1, len(lines)) if lines[i].lstrip().startswith("@C") and "SH2O" in lines[i])
    
    cfmap = parse_header_boundaries(lines[c_header_idx])
    
    sh2o_vals = []
    for i in range(c_header_idx + 1, len(lines)):
        ln = lines[i]
        if ln.startswith("*") or ln.lstrip().startswith("@"):
            break
        if ln.lstrip().startswith("!") or not ln.strip():
            continue
        if ln.strip() and ln.strip()[0].isdigit():
            sh2o_val = float(ln[cfmap["SH2O"][0]:cfmap["SH2O"][1]].strip())
            sh2o_vals.append(sh2o_val)
            
    assert len(sh2o_vals) == 3  # YUKU2101.HMX has 3 layers: 5cm, 100cm, 150cm
    # original value was .464
    # .464 * 1.25 = 0.58
    for val in sh2o_vals:
        assert abs(val - 0.58) < 1e-4


def test_edit_filex_generic_sections_rows_codes_and_initial_conditions(tmp_path):
    filex = tmp_path / "SYNTH.HMX"
    _synthetic_filex(filex)

    edit_filex(
        filex,
        {
            "PPOP": 42.0,
            "row2_spacing": {"field": "PLRS", "value": 44.0, "row": 2, "required": True},
            "irrig_code": {
                "section": "IRRIGATION",
                "field": "IROP",
                "value": "IR005",
                "type": "code",
                "row": 1,
                "required": True,
            },
            "irrig_amount": {
                "section": "IRRIGATION",
                "field": "IRVAL",
                "value": 2.0,
                "op": "mult",
                "row": 2,
                "required": True,
            },
            "fert_n": {
                "section": "FERTILIZERS",
                "header_prefix": "@F",
                "field": "FAMN",
                "value": 45.0,
                "required": True,
            },
        },
        {
            "initial_soil_water_mult": {"value": 1.5, "required": True},
            "soil_nh4": {"field": "SNH4", "value": 0.08, "row": 1, "required": True},
            "soil_no3_add": {"field": "SNO3", "value": 0.10, "op": "add", "required": True},
        },
    )

    lines = filex.read_text(encoding="utf-8").splitlines()

    assert float(_cell(lines, "*PLANTING DETAILS", "@P", "PPOP", 1)) == 42.0
    assert float(_cell(lines, "*PLANTING DETAILS", "@P", "PPOP", 2)) == 42.0
    assert float(_cell(lines, "*PLANTING DETAILS", "@P", "PLRS", 1)) == 50.0
    assert float(_cell(lines, "*PLANTING DETAILS", "@P", "PLRS", 2)) == 44.0

    assert float(_cell(lines, "*INITIAL CONDITIONS", "@C", "SH2O", 1)) == pytest.approx(0.3)
    assert float(_cell(lines, "*INITIAL CONDITIONS", "@C", "SH2O", 2)) == pytest.approx(0.375)
    assert float(_cell(lines, "*INITIAL CONDITIONS", "@C", "SNH4", 1)) == pytest.approx(0.08)
    assert float(_cell(lines, "*INITIAL CONDITIONS", "@C", "SNH4", 2)) == pytest.approx(0.02)
    assert float(_cell(lines, "*INITIAL CONDITIONS", "@C", "SNO3", 1)) == pytest.approx(1.3)
    assert float(_cell(lines, "*INITIAL CONDITIONS", "@C", "SNO3", 2)) == pytest.approx(1.5)

    assert _cell(lines, "*IRRIGATION", "@I", "IROP", 1) == "IR005"
    assert _cell(lines, "*IRRIGATION", "@I", "IROP", 2) == "IR001"
    assert float(_cell(lines, "*IRRIGATION", "@I", "IRVAL", 1)) == 10.0
    assert float(_cell(lines, "*IRRIGATION", "@I", "IRVAL", 2)) == 30.0
    assert float(_cell(lines, "*FERTILIZERS", "@F", "FAMN", 1)) == 45.0


def test_edit_filex_fields_station_and_soil_codes(tmp_path):
    filex = tmp_path / "CNKU2101.HMX"
    filex.write_text(
        "\n".join([
            "*FIELDS",
            "@L ID_FIELD WSTA....  FLSA  ID_SOIL    FLNAME",
            " 1 CTRA2101 CTRA2101   -99  IBSB910015 -99",
        ]) + "\n",
        encoding="utf-8",
    )

    edit_filex(
        filex,
        {
            "id_field": {
                "section": "FIELDS",
                "field": "ID_FIELD",
                "value": "CNKU2101",
                "type": "code",
                "required": True,
            },
            "wsta": {
                "section": "FIELDS",
                "field": "WSTA",
                "value": "CNKU2101",
                "type": "code",
                "required": True,
            },
            "soil": {
                "section": "FIELDS",
                "field": "ID_SOIL",
                "value": "YUKU2101",
                "type": "code",
                "required": True,
            },
        },
        {},
    )

    lines = filex.read_text(encoding="utf-8").splitlines()
    assert _cell(lines, "*FIELDS", "@L", "ID_FIELD", 1) == "CNKU2101"
    assert _cell(lines, "*FIELDS", "@L", "WSTA", 1) == "CNKU2101"
    assert re.search(r"YUKU2101\s+-99", lines[2])
    assert "YUKU2101015" not in lines[2]
    assert "IBSB910015" not in lines[2]


def test_edit_filex_required_generic_field_fails(tmp_path):
    filex = tmp_path / "SYNTH.HMX"
    _synthetic_filex(filex)

    with pytest.raises(KeyError, match="NOT_A_FIELD"):
        edit_filex(
            filex,
            {"bad": {"section": "IRRIGATION", "field": "NOT_A_FIELD", "value": 1, "required": True}},
            {},
        )
