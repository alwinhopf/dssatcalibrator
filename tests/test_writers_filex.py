"""Tests for the FileX management and initial conditions writers."""
import shutil
from pathlib import Path

from dssatcalibrator.writers import edit_filex, parse_header_boundaries

REPO = Path(__file__).resolve().parents[1]
SMOKE_HMX = REPO / "_smoke" / "YUKU2101.HMX"


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
