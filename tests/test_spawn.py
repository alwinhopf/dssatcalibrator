"""End-to-end spawn tests: a real DSSAT hemp run driven by the framework."""
from pathlib import Path

import pytest

from dssatcalibrator.config import load_config, crop_for, active_parameters, resolve_exe
from dssatcalibrator.spawn import (
    _filex_overrides_for,
    _partition_theta,
    parse_cultivars,
    parse_treatments,
    spawn_and_run,
    theta_hash,
)

REPO = Path(__file__).resolve().parents[1]
HEMP_CFG = REPO / "config_hemp.yaml"


@pytest.fixture
def hemp_setup(hemp_dir):
    cfg = load_config(HEMP_CFG)
    # The maintained example config is intentionally machine-neutral. This
    # integration fixture has already discovered a local, optional hemp DSSAT
    # installation, so inject that installation explicitly for the real run.
    dssat_root = hemp_dir.parent
    local_executables = [
        dssat_root / "dscsm048_compiled_4.8.2.with_HM_code.exe",
        dssat_root / "dscsm048",
        dssat_root / "DSCSM048.EXE",
        dssat_root / "dscsm048.exe",
    ]
    local_exe = next((path for path in local_executables if path.is_file()), None)
    if local_exe is None:
        pytest.skip(f"DSSAT hemp executable not found under: {dssat_root}")
    cfg["calibrator"]["dssat_dir"] = str(dssat_root)
    cfg["calibrator"]["dssat_exe"] = str(local_exe)
    cfg["source"]["hemp_dir"] = str(hemp_dir)
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


def test_parse_cultivars(hemp_dir):
    assert parse_cultivars(hemp_dir / "YUKU2101.HMX") == ["IB0008"]
    assert parse_cultivars(hemp_dir / "CNKU2101.HMX") == ["IB0002"]


def test_partition_theta_routes_cultivar_scoped_values():
    specs = [
        {
            "name": "EM-FL__IB0008",
            "base_name": "EM-FL",
            "group": "genetic_cultivar",
            "scope": "cultivar",
            "cultivar": "IB0008",
        },
        {
            "name": "EM-FL__IB0002",
            "base_name": "EM-FL",
            "group": "genetic_cultivar",
            "scope": "cultivar",
            "cultivar": "IB0002",
        },
    ]
    theta = {"EM-FL__IB0008": 28.0, "EM-FL__IB0002": 18.0}
    groups = _partition_theta(theta, specs, exp_id="YUKU2101", cultivars=["IB0008"])
    assert groups == {"genetic_cultivar_by_cultivar": {"IB0008": {"EM-FL": 28.0}}}


def test_theta_hash_supports_filex_code_values():
    value = theta_hash({"irrig_code": "IR004", "x": 1.0})
    assert len(value) == 16
    assert value == theta_hash({"x": 1.0, "irrig_code": "IR004"})
    assert value != theta_hash({"irrig_code": "IR004", "x": 1.0000001})
    assert theta_hash({"x": 1, "irrig_code": "IR004"}) == "06a431c288780c62"
    assert theta_hash({"irrig_code": "IR005", "x": 1.0}) != "31684502a7"


def test_filex_overrides_for_experiment():
    cfg = {
        "filex_overrides": {
            "all": [{"section": "FIELDS", "field": "ID_SOIL", "value": "BASE"}],
            "CNKU2101": [{"section": "FIELDS", "field": "WSTA", "value": "CNKU2101"}],
        }
    }

    assert _filex_overrides_for(cfg, "CNKU2101") == [
        {"section": "FIELDS", "field": "ID_SOIL", "value": "BASE"},
        {"section": "FIELDS", "field": "WSTA", "value": "CNKU2101"},
    ]
    assert _filex_overrides_for(cfg, "YUKU2101") == [
        {"section": "FIELDS", "field": "ID_SOIL", "value": "BASE"},
    ]


def test_spawn_default_reproduces_smoke(hemp_setup, tmp_path):
    cfg, exe, crop, specs, start = hemp_setup
    res = spawn_and_run(start, exp_id="YUKU2101", cfg=cfg, crop=crop,
                        param_specs=specs, run_root=tmp_path, treatments=[1], exe=Path(exe))
    assert res.status == "success", res.message
    adap = res.evaluate[(res.evaluate.variable == "ADAP") & (res.evaluate.treatment == 1)].iloc[0]
    # default coefficients reproduce the installed genotype without a formatting-induced shift
    assert abs(adap["sim"] - 79) <= 1
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
