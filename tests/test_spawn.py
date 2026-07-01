"""End-to-end spawn tests: a real DSSAT hemp run driven by the framework."""
from pathlib import Path

import pytest

from dssatcalibrator.config import load_config, crop_for, active_parameters, resolve_exe
from dssatcalibrator.spawn import spawn_and_run, parse_treatments, theta_hash

REPO = Path(__file__).resolve().parents[1]
HEMP_CFG = REPO / "config_hemp.yaml"


@pytest.fixture
def hemp_setup(hemp_dir):
    cfg = load_config(HEMP_CFG)
    exe = resolve_exe(cfg)
    if not Path(exe).exists():
        pytest.skip(f"DSSAT exe not found: {exe}")
    crop = crop_for(cfg, "HM")
    specs = active_parameters(cfg)
    start = {p["name"]: p["start"] for p in specs}
    return cfg, exe, crop, specs, start


def test_parse_treatments(hemp_dir):
    trts = parse_treatments(hemp_dir / "YUKU2101.HMX")
    assert trts == [1, 2, 3, 4]


def test_theta_hash_supports_filex_code_values():
    assert theta_hash({"irrig_code": "IR004", "x": 1.0}) == "31684502a7"
    assert theta_hash({"x": 1, "irrig_code": "IR004"}) == "31684502a7"
    assert theta_hash({"irrig_code": "IR005", "x": 1.0}) != "31684502a7"


def test_spawn_default_reproduces_smoke(hemp_setup, tmp_path):
    cfg, exe, crop, specs, start = hemp_setup
    res = spawn_and_run(start, exp_id="YUKU2101", cfg=cfg, crop=crop,
                        param_specs=specs, run_root=tmp_path, treatments=[1], exe=Path(exe))
    assert res.status == "success", res.message
    adap = res.evaluate[(res.evaluate.variable == "ADAP") & (res.evaluate.treatment == 1)].iloc[0]
    # default coefficients reproduce the smoke run (anthesis sim 76 vs measured 75)
    assert abs(adap["sim"] - 76) <= 1
    assert adap["meas"] == 75
    # time-series parsed with biomass present
    assert (res.plantgro["CWAD"] > 0).any()


def test_spawn_shifts_phenology(hemp_setup, tmp_path):
    cfg, exe, crop, specs, start = hemp_setup
    base = spawn_and_run(start, exp_id="YUKU2101", cfg=cfg, crop=crop,
                         param_specs=specs, run_root=tmp_path / "base",
                         treatments=[1], exe=Path(exe))
    early = dict(start)
    early["EM-FL"] = 14.0   # much shorter emergence->flower => earlier anthesis
    shifted = spawn_and_run(early, exp_id="YUKU2101", cfg=cfg, crop=crop,
                            param_specs=specs, run_root=tmp_path / "early",
                            treatments=[1], exe=Path(exe))
    a0 = base.evaluate.query("variable=='ADAP' and treatment==1").iloc[0]["sim"]
    a1 = shifted.evaluate.query("variable=='ADAP' and treatment==1").iloc[0]["sim"]
    assert a1 < a0, f"expected earlier anthesis, got {a1} vs {a0}"
