"""Tests for config loading and active-parameter enumeration."""
from pathlib import Path

from dssatcalibrator import config as cfgmod

REPO = Path(__file__).resolve().parents[1]
HEMP_CFG = REPO / "config_hemp.yaml"


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
