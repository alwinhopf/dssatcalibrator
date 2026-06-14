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
        "treatment": [1, 1], "date": DATES, "CWAD": cwad_sim,
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
