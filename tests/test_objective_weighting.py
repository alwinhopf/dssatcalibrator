"""Tests for the new objective behaviour: count_scale vs unified, and obs_autocorr."""
from types import SimpleNamespace

import numpy as np
import pandas as pd

from dssatcalibrator import objective as obj
from dssatcalibrator.observations import SCHEMA

DATES = [pd.Timestamp("2021-07-01") + pd.Timedelta(days=i) for i in range(5)]

CFG = {
    "engine": {"timeseries_outputs": {"biomass": "CWAD"},
               "scalar_outputs": {"grain_yield": "HWAM"}},
    "objective": {"weighting": "unified", "weights": {"biomass": 1.0, "grain_yield": 1.0},
                  "error_model": {"biomass": {"type": "absolute", "value": 100.0},
                                  "grain_yield": {"type": "absolute", "value": 100.0}}},
}


def _results(cwad_sim):
    ev = pd.DataFrame({"treatment": [1], "variable": ["HWAM"], "sim": [1000.0], "meas": [1000.0]})
    pg = pd.DataFrame({"treatment": [1] * 5, "date": DATES, "CWAD": cwad_sim})
    return {"E1": SimpleNamespace(evaluate=ev, plantgro=pg)}


def _obs():
    vals = [100.0, 200.0, 300.0, 400.0, 500.0]   # monotonic -> strongly autocorrelated
    rows = [("E1", 1, "CWAD", "timeseries", d, v, np.nan, 1.0) for d, v in zip(DATES, vals)]
    return pd.DataFrame(rows, columns=SCHEMA)


def test_count_scale_differs_from_unified():
    # biomass has 5 points, yield has 1 -> averaging (count_scale) vs summing
    # (unified) over the two variable groups must give different totals.
    results = _results([150.0, 250.0, 350.0, 450.0, 550.0])   # biomass off by 50 each, yield exact
    u = obj.score(results, _obs(), {**CFG, "objective": {**CFG["objective"], "weighting": "unified"}})
    c = obj.score(results, _obs(), {**CFG, "objective": {**CFG["objective"], "weighting": "count_scale"}})
    assert np.isfinite(u.score) and np.isfinite(c.score)
    assert u.score != c.score


def test_obs_autocorr_downweights_timeseries():
    simulated = [110.0, 220.0, 330.0, 440.0, 550.0]
    base = obj.build_residuals(_results(simulated), _obs(), CFG)
    cfg2 = {**CFG, "objective": {**CFG["objective"], "obs_autocorr": True}}
    aut = obj.build_residuals(_results(simulated), _obs(), cfg2)
    w_base = base[base["kind"] == "timeseries"]["weight"].sum()
    w_aut = aut[aut["kind"] == "timeseries"]["weight"].sum()
    assert w_aut < w_base        # serial-correlation down-weighting shrank the series


def test_row_level_weights_affect_unified_score():
    obs = _obs()
    obs.loc[0, "weight"] = 0.01
    sim = [300.0, 200.0, 300.0, 400.0, 500.0]

    downweighted = obj.score(_results(sim), obs, CFG)
    equal = obj.score(_results(sim), _obs(), CFG)

    assert downweighted.score < equal.score
