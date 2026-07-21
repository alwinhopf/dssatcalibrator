"""Offline tests for the in-season forecast + new-crop scaffolding features.

Covers the pure-Python logic; anything needing DSSAT runs / live network is not
exercised here (those are integration-tier).
"""

from datetime import date
import os
from pathlib import Path
import stat
from types import SimpleNamespace
import sys

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
SMOKE = REPO / "_smoke"
if not (SMOKE / "HMGRO048.CUL").exists():
    SMOKE = REPO.parent / "DSSAT48" / "Genotype"


# --------------------------------------------------------------------------- #
# weather.py
# --------------------------------------------------------------------------- #
def test_write_and_fill_weather(tmp_path):
    from dssatcalibrator.weather import build_provider, read_wth, write_wth, fill_gap, WTH_COLS

    df = pd.DataFrame({
        "date": pd.date_range("2024-05-01", periods=10, freq="D"),
        "SRAD": np.linspace(18, 22, 10), "TMAX": np.linspace(28, 32, 10),
        "TMIN": np.linspace(15, 18, 10), "RAIN": np.zeros(10),
    })
    p = write_wth(tmp_path / "TEST.WTH", "TEST", 14.0, -91.0, df)
    text = p.read_text()
    assert "@DATE  SRAD  TMAX  TMIN  RAIN" in text
    assert "24122" in text  # 2024-05-01 is DOY 122
    reread = read_wth(p)
    assert list(reread.columns) == ["date", *WTH_COLS]
    assert reread.loc[0, "date"] == pd.Timestamp("2024-05-01")
    assert reread.loc[0, "SRAD"] == pytest.approx(18.0)
    assert build_provider({"weather": {"provider": "dssatutils"}}).__class__.__name__ == "DssatutilsWeatherProvider"

    filled = fill_gap(df, "2024-05-20", method="persistence")
    assert filled["filled"].sum() == 10          # 10 appended days
    assert filled["date"].max() == pd.Timestamp("2024-05-20")
    clim = fill_gap(df, "2024-05-15", method="climatology")
    assert clim["filled"].any()
    assert fill_gap(df, "2024-05-20", method="none")["filled"].sum() == 0


def test_nasa_power_parse():
    from dssatcalibrator.weather import NasaPowerProvider
    payload = {"properties": {"parameter": {
        "ALLSKY_SFC_SW_DWN": {"20240501": 20.5, "20240502": 21.0},
        "T2M_MAX": {"20240501": 30.0, "20240502": 31.0},
        "T2M_MIN": {"20240501": 18.0, "20240502": -99},
        "PRECTOTCORR": {"20240501": 0.0, "20240502": 5.0},
    }}}
    df = NasaPowerProvider._parse_power(payload)
    assert list(df["date"]) == [pd.Timestamp("2024-05-01"), pd.Timestamp("2024-05-02")]
    assert df.loc[0, "SRAD"] == 20.5
    assert np.isnan(df.loc[1, "TMIN"])           # -99 -> NaN


# --------------------------------------------------------------------------- #
# forecast.py
# --------------------------------------------------------------------------- #
def test_ensemble_percentiles_and_anchor():
    from dssatcalibrator.forecast import ensemble_percentiles, anchor_correction, lead_time_table

    dates = pd.date_range("2024-06-01", periods=5, freq="D")
    curves = [pd.DataFrame({"date": dates, "LAID": np.full(5, v)}) for v in (1.0, 2.0, 3.0)]
    pct = ensemble_percentiles(curves, "LAID")
    assert np.isclose(pct.loc[0, "p50"], 2.0)
    assert np.isclose(pct.loc[0, "mean"], 2.0)
    assert (pct["n"] == 3).all()

    # anchor: obs=2.5 at day0 where p50=2.0 -> +0.5 correction decaying over 4 days
    adj = anchor_correction(pct, last_obs_value=2.5, last_obs_date=dates[0], decay_days=4)
    assert np.isclose(adj.loc[0, "p50_adj"], 2.5)            # fully corrected at anchor
    assert np.isclose(adj.loc[4, "p50_adj"], 2.0)            # decayed back to model
    assert adj.loc[0, "anchor_weight"] == 1.0

    lt = lead_time_table(pct, dates[0])
    assert lt.loc[0, "lead_days"] == 0
    # p90-p10 of [1,2,3] with linear interpolation = 2.8 - 1.2 = 1.6
    assert np.isclose(lt.loc[0, "spread"], 1.6)


def test_behavioural_thetas():
    from dssatcalibrator.forecast import _behavioural_thetas
    space = SimpleNamespace(names=["A", "B"])
    design = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6], "score": [0.3, 0.1, 0.2]})
    result = SimpleNamespace(space=space, design=design, best_theta={"A": 2, "B": 5})
    thetas = _behavioural_thetas(result, n=2)
    assert thetas[0] == {"A": 2, "B": 5}                    # best by score first
    # n=0 -> best only
    assert _behavioural_thetas(result, n=0) == [{"A": 2, "B": 5}]


# --------------------------------------------------------------------------- #
# diagnostics.py
# --------------------------------------------------------------------------- #
def _fake_result():
    cfg = {"parameters": {"g": {
        "A": {"active": True, "min": 0, "max": 10, "start": 5},
        "B": {"active": True, "min": 0, "max": 10, "start": 5, "prior": {"dist": "normal", "sd": 2}},
    }}}
    space = SimpleNamespace(names=["A", "B"])
    rng = np.random.default_rng(0)
    # A is well constrained (tight), B is wide (~prior)
    design = pd.DataFrame({
        "A": rng.normal(5, 0.2, 50), "B": rng.uniform(0, 10, 50),
        "score": rng.uniform(0, 1, 50),
    })
    best = SimpleNamespace(per_var={"LAI": {"EF": 0.8, "nRMSE_pct": 12, "n": 9},
                                    "yield": {"EF": -0.5, "nRMSE_pct": 70, "n": 3}})
    return SimpleNamespace(cfg=cfg, space=space, design=design, best=best,
                           best_theta={"A": 5, "B": 5})


def test_identifiability_and_structural():
    from dssatcalibrator.diagnostics import identifiability, structural_adequacy
    ident = identifiability(_fake_result())
    a = ident.set_index("parameter").loc["A"]
    b = ident.set_index("parameter").loc["B"]
    assert a["identifiable"] and not b["identifiable"]      # A tight, B ~ prior

    struct = structural_adequacy(_fake_result()).set_index("variable")
    assert not struct.loc["LAI", "flag"]
    assert struct.loc["yield", "flag"]                      # EF<0 -> flagged


# --------------------------------------------------------------------------- #
# writers .SPE + scaffold
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not (SMOKE / "HMGRO048.SPE").exists(), reason="needs bundled SPE")
def test_edit_species(tmp_path):
    from dssatcalibrator.writers import edit_species
    spe = tmp_path / "X.SPE"
    spe.write_text((SMOKE / "HMGRO048.SPE").read_text(errors="replace"))
    edit_species(spe, {"PARMAX,PHTMAX,KCAN": 42.0})
    line = next(l for l in spe.read_text().splitlines()
                if "PARMAX,PHTMAX,KCAN" in l and not l.lstrip().startswith("!"))
    assert "42" in line.split()[0]

    with pytest.raises(KeyError):                            # no match
        edit_species(spe, {"NONEXISTENT_KEY_XYZ": 1.0})


@pytest.mark.skipif(not (SMOKE / "HMGRO048.CUL").exists(), reason="needs bundled CUL")
def test_cul_bounds_and_scaffold_yaml():
    from dssatcalibrator.writers import read_cul_calibration_bounds
    from dssatcalibrator.scaffold import _emit_parameters_yaml, _role
    # smoke file may lack 999991/999992 rows -> returns {} without error
    bounds = read_cul_calibration_bounds(SMOKE / "HMGRO048.CUL")
    assert isinstance(bounds, dict)

    coeffs = {"CSDL": {"min": 11.0, "max": 16.0, "start": 12.8, "role": _role("CSDL")},
              "LFMAX": {"min": 0.5, "max": 2.0, "start": 1.4, "role": _role("LFMAX")}}
    y = _emit_parameters_yaml(coeffs)
    assert "genetic_cultivar" in y
    assert "CSDL" in y and "role: obligatory" in y          # phenology -> obligatory
    assert "LFMAX" in y and "role: candidate" in y          # growth -> candidate


# --------------------------------------------------------------------------- #
# satellite cloud masking + obs operator
# --------------------------------------------------------------------------- #
def test_satellite_cloud_mask_and_operator(tmp_path):
    from dssatcalibrator.sources.satellite import SentinelLAISource
    csv = tmp_path / "lai.csv"
    pd.DataFrame({
        "exp_id": ["E1", "E1", "E1"],
        "treatment": [1, 1, 1],
        "date": ["2024-06-01", "2024-06-08", "2024-06-15"],
        "value": [2.0, 3.0, 4.0],
        "cloud_fraction": [0.1, 0.9, 0.0],
    }).to_csv(csv, index=False)

    src = SentinelLAISource({"data_path": str(csv), "max_cloud_fraction": 0.5,
                             "obs_operator": {"scale": 0.9, "offset": 0.1}})
    df = src.fetch("E1", (date(2024, 5, 1), date(2024, 7, 1)))
    assert len(df) == 2                                     # cloudy 0.9 dropped
    # operator applied: 0.9*2.0 + 0.1 = 1.9
    assert np.isclose(df.iloc[0]["value"], 1.9)


# --------------------------------------------------------------------------- #
# planting dates + orchestrator helpers + config defaults
# --------------------------------------------------------------------------- #
def test_planting_dates_from_obs():
    from dssatcalibrator.observations import Observations
    tbl = pd.DataFrame({
        "exp_id": ["E1", "E2", "E1"],
        "variable": ["planting_date", "sowing_date", "LAID"],
        "date": [pd.Timestamp("2024-05-01"), pd.Timestamp("2024-05-10"), pd.Timestamp("2024-06-01")],
        "value": [np.nan, np.nan, 2.0],
    })
    pd_ = Observations(tbl).planting_dates()
    assert pd_["E1"] == pd.Timestamp("2024-05-01")
    assert pd_["E2"] == pd.Timestamp("2024-05-10")


def test_folds_and_staging():
    from dssatcalibrator.orchestrator import _make_folds, _apply_staging, _year_key, _site_key
    exps = ["YUKU2101", "YUKU2201", "YUFE2101"]
    assert _year_key("YUKU2101") == "21"
    assert _site_key("YUKU2101") == "YUKU"
    loeo = _make_folds(exps, "loeo", 42)
    assert len(loeo) == 3
    years = dict(_make_folds(exps, "year", 42))
    assert set(years) == {"year_21", "year_22"}             # YUKU2101/YUFE2101 vs YUKU2201

    cfg = {"method": {"staging": {"freeze_groups": ["genetic_ecotype"], "freeze_params": ["XFRT"]}},
           "parameters": {
               "genetic_cultivar": {"CSDL": {"active": True}, "XFRT": {"active": True}},
               "genetic_ecotype": {"PL-EM": {"active": True}}}}
    out = _apply_staging(cfg)
    assert out["parameters"]["genetic_cultivar"]["CSDL"]["active"] is True
    assert out["parameters"]["genetic_cultivar"]["XFRT"]["active"] is False
    assert out["parameters"]["genetic_ecotype"]["PL-EM"]["active"] is False


def test_config_defaults_present():
    from dssatcalibrator.config import DEFAULTS
    assert DEFAULTS["weather"]["provider"] == "file"
    assert DEFAULTS["execution"]["backend"] == "native"
    assert DEFAULTS["soil"]["provider"] == "file"
    assert DEFAULTS["forecast"]["active"] is False
    assert DEFAULTS["gating"]["species"] == "blocked"
    assert DEFAULTS["management_options"]["use_source_planting_date"] is False


def test_shared_execution_backend_uses_dssatengine(monkeypatch, tmp_path):
    from dssatcalibrator import spawn as sp

    calls = {}

    def fake_normalize(start, end, treatment_list=None, treatments=None):
        calls["normalize"] = list(treatment_list)
        return [5, 1]

    def fake_write(filex, trts, batch_path, run_mode="experiment"):
        calls["write"] = (filex, list(trts), run_mode)
        Path(batch_path).write_text("batch", encoding="utf-8")

    def fake_run(run_dir, exe, run_mode_flag="A", filex="", model=None, timeout=None):
        calls["run"] = (run_dir, exe, run_mode_flag, model, timeout)

    fake_engine = SimpleNamespace(
        normalize_treatment_list=fake_normalize,
        write_dssbatch=fake_write,
        run_dssat=fake_run,
    )
    monkeypatch.setitem(sys.modules, "dssatengine", fake_engine)

    treatments = sp._normalize_treatments([5, "1", 5], "dssatengine")
    assert treatments == [5, 1]
    sp._write_batch(tmp_path, "TEST.HMX", treatments, "dssatengine")
    err = sp._run_backend_dssat(tmp_path, Path("fake_exe"), {"model": "CRGRO048"},
                                "dssatengine", timeout=12)

    assert err == ""
    assert calls["normalize"] == [5, "1", 5]
    assert calls["write"] == ("TEST.HMX", [5, 1], "experiment")
    assert calls["run"][2:] == ("B", "CRGRO048", 12)


def _fake_dssat_exe(tmp_path: Path, exit_code: int = 0) -> Path:
    if os.name == "nt":
        exe = tmp_path / "fake_dssat.bat"
        exe.write_text(
            "@echo off\n"
            "echo %* > native_args.txt\n"
            "echo native dssat stdout\n"
            f"exit /B {exit_code}\n",
            encoding="utf-8",
        )
    else:
        exe = tmp_path / "fake_dssat"
        exe.write_text(
            "#!/bin/sh\n"
            "echo \"$@\" > native_args.txt\n"
            "echo native dssat stdout\n"
            f"exit {exit_code}\n",
            encoding="utf-8",
        )
        exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    return exe


def test_native_execution_backend_captures_logs_and_uses_uppercase_batch(tmp_path):
    from dssatcalibrator import spawn as sp

    exe = _fake_dssat_exe(tmp_path)
    err = sp._run_backend_dssat(tmp_path, exe, {"model": "CRGRO048"},
                                "native", timeout=10)

    assert err == ""
    assert (tmp_path / "native_args.txt").read_text(encoding="utf-8").strip() == (
        "CRGRO048 B DSSBatch.V48"
    )
    assert "native dssat stdout" in (
        tmp_path / "dssat_B_stdout_stderr.log"
    ).read_text(encoding="utf-8")


def test_native_execution_backend_reports_nonzero_status(tmp_path):
    from dssatcalibrator import spawn as sp

    exe = _fake_dssat_exe(tmp_path, exit_code=9)
    err = sp._run_backend_dssat(tmp_path, exe, {"model": "CRGRO048"},
                                "native", timeout=10)

    assert "status 9" in err
    assert "native dssat stdout" in err
