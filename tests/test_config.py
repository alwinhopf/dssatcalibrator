"""Tests for config loading and active-parameter enumeration."""
from pathlib import Path

import pytest

from dssatcalibrator import config as cfgmod

REPO = Path(__file__).resolve().parents[1]
HEMP_CFG = REPO / "config_hemp.yaml"


def test_validate_config_accepts_reference_config():
    # The shipped reference config must always pass validation.
    cfg = cfgmod.load_config(HEMP_CFG)            # load_config validates by default
    assert cfgmod.validate_config(cfg) is cfg


def test_validate_config_collects_all_problems():
    bad = {
        "method": {"preset": "Z", "validation": {"scheme": "weekly"}},
        "objective": {"weighting": "bogus"},
        "execution": {"backend": "spark"},
        "assimilation": {"mode": "kalman"},
        "gating": {"cultivar": "loose"},
        "calibrator": {"num_cores": -2},
        "parameters": {
            "g": {
                "A": {"active": True, "min": 5, "max": 1},                # inverted bounds
                "B": {"active": True, "min": 0, "max": 10, "start": 50},  # start out of bounds
                "C": {"active": True, "min": 0, "max": 1, "prior": {"dist": "cauchy"}},
            }
        },
    }
    with pytest.raises(ValueError) as ei:
        cfgmod.validate_config(bad)
    msg = str(ei.value)
    for needle in ("method.preset", "validation.scheme", "objective.weighting",
                   "execution.backend", "assimilation.mode", "gating.cultivar",
                   "num_cores", "min (5) must be < max (1)", "start (50) is outside",
                   "prior.dist 'cauchy'"):
        assert needle in msg, needle


def test_validate_config_requires_an_active_parameter():
    cfg = {"parameters": {"g": {"A": {"active": False, "min": 0, "max": 1}}}}
    with pytest.raises(ValueError, match="No active parameters"):
        cfgmod.validate_config(cfg)


def test_load_config_validate_false_skips_check(tmp_path):
    # A config with no active parameters loads fine when validation is disabled.
    p = tmp_path / "c.yaml"
    p.write_text("parameters: {}\n")
    cfg = cfgmod.load_config(p, validate=False)
    assert cfg["parameters"] == {}


def test_load_hemp_config():
    cfg = cfgmod.load_config(HEMP_CFG)
    assert cfg["calibrator"]["name"] == "hemp_yunnan_2021_2022"
    assert cfg["method"]["preset"] == "A"
    crop = cfgmod.crop_for(cfg, "HM")
    assert crop["model"] == "CRGRO048"
    assert crop["cultivar_anchor"] == "IB0008"


def test_active_parameters():
    cfg = cfgmod.load_config(HEMP_CFG)
    active = cfgmod.active_parameters(cfg)
    names = {p["name"] for p in active}
    # updated with additional cultivar and ecotype parameters
    assert names == {"CSDL", "PPSEN", "EM-FL", "FL-SD", "SD-PM", "LFMAX", "SLAVR", "SIZLF", "WTPSD", "SFDUR", "SDPDV", "PL-EM", "RHGHT"}
    for p in active:
        assert p["min"] < p["max"]
        assert p["min"] <= p["start"] <= p["max"], p["name"]


def test_env_override(monkeypatch):
    monkeypatch.setenv("DSSATCAL_NUM_CORES", "5")
    cfg = cfgmod.load_config(HEMP_CFG)
    assert cfg["calibrator"]["num_cores"] == 5


def test_resolve_exe():
    cfg = cfgmod.load_config(HEMP_CFG)
    exe = cfgmod.resolve_exe(cfg)
    assert "with_HM_code" in str(exe)


def test_shared_stack_defaults_and_template_env(monkeypatch, tmp_path):
    assert cfgmod.DEFAULTS["execution"]["backend"] == "native"
    assert cfgmod.DEFAULTS["soil"]["provider"] == "file"

    monkeypatch.setenv("DSSAT_TEMPLATE_DIR", str(tmp_path))
    assert cfgmod.resolve_template_dir({}) == tmp_path

    cfg = cfgmod.load_config(HEMP_CFG)
    paths = cfgmod.resolve_dssat_paths(cfg)
    assert paths["genotype"].name == "Genotype"
    assert paths["weather"].name == "Weather"
    assert paths["soil"].name == "Soil"
