"""Tests for the objective (alignment, metrics, weighting) using synthetic data."""
from types import SimpleNamespace

import numpy as np
import pandas as pd

from dssatcalibrator import objective as obj
from dssatcalibrator.observations import SCHEMA

CFG = {
    "engine": {
        "timeseries_outputs": {"biomass": "CWAD"},
        "scalar_outputs": {"anthesis": "ADAP", "grain_yield": "HWAM"},
    },
    "objective": {
        "weighting": "unified",
        "weights": {"biomass": 1.0, "anthesis": 1.0, "grain_yield": 1.0},
        "error_model": {
            "anthesis": {"type": "absolute", "value": 3},
            "grain_yield": {"type": "relative", "value": 0.1},
            "biomass": {"type": "relative", "value": 0.2},
        },
    },
}

DATES = [pd.Timestamp("2021-09-17"), pd.Timestamp("2021-10-05")]


def _result(adap_sim, hwam_sim, cwad_sim):
    ev = pd.DataFrame({
        "treatment": [1, 1], "variable": ["ADAP", "HWAM"],
        "sim": [adap_sim, hwam_sim], "meas": [75.0, 1000.0],
    })
    pg = pd.DataFrame({
        "treatment": [1, 1], "date": DATES, "DAP": [75, 93], "CWAD": cwad_sim,
    })
    return SimpleNamespace(evaluate=ev, plantgro=pg)


def _obs():
    rows = [
        ("E1", 1, "CWAD", "timeseries", DATES[0], 5000.0, np.nan, 1.0),
        ("E1", 1, "CWAD", "timeseries", DATES[1], 8000.0, np.nan, 1.0),
    ]
    return pd.DataFrame(rows, columns=SCHEMA)


def test_metrics_perfect():
    m = obj.metrics([1, 2, 3], [1, 2, 3])
    assert m["RMSE"] == 0
    assert m["d"] == 1 and m["EF"] == 1 and round(m["R2"], 6) == 1
    assert m["n"] == 3


def test_build_residuals_both_paths():
    results = {"E1": _result(75.0, 1000.0, [5000.0, 8000.0])}
    resid = obj.build_residuals(results, _obs(), CFG)
    # 2 scalars (ADAP, HWAM) + 2 time-series (CWAD) = 4 rows
    assert len(resid) == 4
    assert set(resid["user_var"]) == {"anthesis", "grain_yield", "biomass"}
    assert "phenology" in set(resid["kind"]) and "timeseries" in set(resid["kind"])


def test_unconfigured_plantgro_columns_are_not_scored():
    """Available PlantGro columns are not automatically calibration targets."""
    obs = pd.concat([
        _obs(),
        pd.DataFrame([
            ("E1", 1, "RWAD", "timeseries", DATES[0], 250.0, np.nan, 1.0),
        ], columns=SCHEMA),
    ], ignore_index=True)
    result = _result(75.0, 1000.0, [5000.0, 8000.0])
    result.plantgro["RWAD"] = [250.0, 300.0]

    resid = obj.build_residuals({"E1": result}, obs, CFG)

    assert "RWAD" not in set(resid["dssat"])
    assert len(resid[resid["kind"] == "timeseries"]) == 2


def test_filea_phenology_date_maps_to_dap_output():
    obs = pd.DataFrame([
        ("E1", 1, "ADAT", "phenology", DATES[0], 21260.0, np.nan, 1.0),
    ], columns=SCHEMA)

    resid = obj.build_residuals({"E1": _result(75.0, 1000.0, [5000.0, 8000.0])}, obs, CFG)

    anthesis = resid[resid["user_var"] == "anthesis"].iloc[0]
    assert anthesis["dssat"] == "ADAP"
    assert anthesis["obs"] == 75.0
    assert anthesis["sim"] == 75.0


def test_late_timeseries_obs_uses_last_simulated_value():
    ev = pd.DataFrame({
        "treatment": [1, 1], "variable": ["ADAP", "HWAM"],
        "sim": [75.0, 1000.0], "meas": [75.0, 1000.0],
    })
    pg = pd.DataFrame({
        "treatment": [1], "date": [DATES[0]], "DAP": [75], "CWAD": [5000.0],
    })
    results = {"E1": SimpleNamespace(evaluate=ev, plantgro=pg)}

    resid = obj.build_residuals(results, _obs(), CFG)

    biomass = resid[resid["user_var"] == "biomass"].sort_values("date")
    assert len(biomass) == 2
    assert biomass.iloc[-1]["date"] == DATES[1]
    assert biomass.iloc[-1]["sim"] == 5000.0
    assert biomass.iloc[-1]["obs"] == 8000.0


def test_configured_zero_observations_are_ignored():
    obs = pd.DataFrame([
        ("E1", 1, "CWAD", "timeseries", DATES[0], 0.0, np.nan, 1.0),
        ("E1", 1, "CWAD", "timeseries", DATES[1], 8000.0, np.nan, 1.0),
    ], columns=SCHEMA)
    cfg = {
        **CFG,
        "objective": {
            **CFG["objective"],
            "ignore_zero_observations": ["biomass"],
        },
    }

    resid = obj.build_residuals({"E1": _result(75.0, 1000.0, [5000.0, 8000.0])}, obs, cfg)

    biomass = resid[resid["user_var"] == "biomass"]
    assert len(biomass) == 1
    assert biomass.iloc[0]["obs"] == 8000.0


def test_observation_sigma_and_weight_override_objective_defaults():
    obs = _obs().iloc[[0]].copy()
    obs.loc[:, "sigma"] = 50.0
    obs.loc[:, "weight"] = 0.25
    cfg = {
        **CFG,
        "objective": {
            **CFG["objective"],
            "weights": {**CFG["objective"]["weights"], "biomass": 0.5},
        },
    }

    resid = obj.build_residuals(
        {"E1": _result(75.0, 1000.0, [5000.0, 8000.0])}, obs, cfg
    )
    biomass = resid[resid["user_var"] == "biomass"].iloc[0]

    assert biomass["sigma"] == 50.0
    assert biomass["weight"] == 0.25


def test_perfect_fit_scores_zero():
    results = {"E1": _result(75.0, 1000.0, [5000.0, 8000.0])}
    r = obj.score(results, _obs(), CFG)
    assert r.score == 0.0
    assert r.loglik == 0.0
    assert r.per_var["biomass"]["RMSE"] == 0


def test_biased_fit_scores_worse():
    good = obj.score({"E1": _result(75.0, 1000.0, [5000.0, 8000.0])}, _obs(), CFG)
    bad = obj.score({"E1": _result(90.0, 700.0, [3000.0, 6000.0])}, _obs(), CFG)
    assert bad.score > good.score
    assert bad.loglik < good.loglik
    # anthesis off by 15 days with sigma=3 -> large standardized residual
    assert bad.per_var["anthesis"]["MBE"] == 15.0


def test_weighting_modes_run():
    results = {"E1": _result(80.0, 900.0, [4000.0, 7000.0])}
    for mode in ("unified", "sigma", "count_scale", "user"):
        cfg = {**CFG, "objective": {**CFG["objective"], "weighting": mode}}
        r = obj.score(results, _obs(), cfg)
        assert np.isfinite(r.score) and r.score > 0
