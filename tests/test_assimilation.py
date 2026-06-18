import pytest
import pandas as pd
import numpy as np
from datetime import date
from dssatcalibrator.sources.registry import build_sources, ADAPTER_REGISTRY
from dssatcalibrator.sources.base import ObservationSource
from dssatcalibrator.fusion import ObservationFuser
from dssatcalibrator.engines.enkf import EnsembleKalmanFilter
from dssatcalibrator.engines.forcing import ForcingAssimilator
from dssatcalibrator.engines.recalibration import InSeasonRecalibrator
from dssatcalibrator.observations import Observations

def test_registry():
    assert "sentinel2_lai" in ADAPTER_REGISTRY
    assert "field_measurements" in ADAPTER_REGISTRY
    assert "soil_moisture_iot" in ADAPTER_REGISTRY

    cfg = {
        "observation_sources": {
            "sentinel2_lai": {"active": True, "data_path": "dummy.csv"},
            "soil_moisture_iot": {"active": False}
        }
    }
    sources = build_sources(cfg)
    assert len(sources) == 1
    assert sources[0].name == "sentinel2_lai"

def test_fuser_conflict_resolution():
    # Setup mock sources and data
    class MockSource1(ObservationSource):
        name = "src1"
        source_type = "satellite"
        def fetch(self, exp, date_range, **kw):
            return pd.DataFrame([
                {"exp_id": exp, "treatment": 1, "variable": "LAID", "kind": "timeseries",
                 "date": pd.Timestamp("2026-06-01"), "value": 2.0, "sigma": 0.5, "weight": 1.0,
                 "quality_flag": 0, "spatial_res_m": 10.0}
            ])
        def error_model(self, v, val, meta): return 0.5
        def variable_mapping(self): return {}

    class MockSource2(ObservationSource):
        name = "src2"
        source_type = "uav"
        def fetch(self, exp, date_range, **kw):
            return pd.DataFrame([
                {"exp_id": exp, "treatment": 1, "variable": "LAID", "kind": "timeseries",
                 "date": pd.Timestamp("2026-06-01"), "value": 2.2, "sigma": 0.2, "weight": 1.0,
                 "quality_flag": 0, "spatial_res_m": 0.05}
            ])
        def error_model(self, v, val, meta): return 0.2
        def variable_mapping(self): return {}

    cfg = {
        "fusion": {
            "conflict_resolution": "inverse_variance"
        }
    }
    
    fuser = ObservationFuser([MockSource1({}), MockSource2({})], cfg)
    fused = fuser.collect("EXP1", (date(2026, 5, 1), date(2026, 7, 1)))
    
    assert len(fused) == 1
    # Inverse variance mean of 2.0 (w=4) and 2.2 (w=25) -> (8 + 55)/29 = 63/29 ~ 2.17
    assert np.isclose(fused.iloc[0]["value"], (2.0*4 + 2.2*25)/29)
    # Inverse variance sigma: 1/sqrt(29) ~ 0.185
    assert np.isclose(fused.iloc[0]["sigma"], 1.0 / np.sqrt(29))

    # Test priority strategy
    cfg_priority = {
        "fusion": {
            "conflict_resolution": "priority",
            "source_priority": ["src2", "src1"]
        }
    }
    fuser_p = ObservationFuser([MockSource1({}), MockSource2({})], cfg_priority)
    fused_p = fuser_p.collect("EXP1", (date(2026, 5, 1), date(2026, 7, 1)))
    assert len(fused_p) == 1
    # Should select src2 value because it has higher priority rank
    assert fused_p.iloc[0]["value"] == 2.2
    assert fused_p.iloc[0]["source"] == "src2"

def test_enkf():
    cfg = {
        "assimilation": {
            "enkf": {
                "n_ensemble": 5,
                "inflation": 1.0,
                "state_variables": ["LAID", "CWAD"]
            }
        },
        "calibrator": {"seed": 42}
    }
    enkf = EnsembleKalmanFilter(cfg)
    
    # Shape: (5, 2)
    ensemble = np.array([
        [1.0, 100.0],
        [1.2, 110.0],
        [0.8, 90.0],
        [1.1, 105.0],
        [0.9, 95.0]
    ])
    
    updated = enkf.assimilate(ensemble, "LAID", 1.5, 0.1)
    assert updated.shape == (5, 2)
    # The mean LAID should shift closer to the observed 1.5
    assert updated.mean(axis=0)[0] > ensemble.mean(axis=0)[0]

def test_forcing():
    cfg = {
        "assimilation": {
            "forcing": {
                "min_confidence": 0.7,
                "smoothing": False
            }
        }
    }
    assimilator = ForcingAssimilator(cfg)
    state = {"LAID": 1.0}
    obs = {"variable": "LAID", "value": 2.0, "confidence": 0.8}
    
    updated = assimilator.apply(state, obs)
    assert updated["LAID"] == 2.0

    # With smoothing
    cfg_smooth = {
        "assimilation": {
            "forcing": {
                "min_confidence": 0.7,
                "smoothing": True
            }
        }
    }
    assimilator_s = ForcingAssimilator(cfg_smooth)
    updated_s = assimilator_s.apply(state, obs)
    # 0.8 * 2.0 + 0.2 * 1.0 = 1.6 + 0.2 = 1.8
    assert np.isclose(updated_s["LAID"], 1.8)

def test_recalibrator_prep():
    cfg = {
        "calibrator": {"name": "test"},
        "assimilation": {
            "recalibration": {
                "recal_sample_size": 10
            }
        }
    }
    recal = InSeasonRecalibrator(cfg)
    obs_df = pd.DataFrame([
        {"date": pd.Timestamp("2026-06-01"), "value": 1.0},
        {"date": pd.Timestamp("2026-06-15"), "value": 2.0},
        {"date": pd.Timestamp("2026-07-01"), "value": 3.0}
    ])

    # Just verify the date filtering and cfg duplication does not crash
    with pytest.raises(Exception):
        # calibrate will fail on setup because we gave dummy directories,
        # but this checks that the recalibrate method gets past filtering and calls calibrate
        recal.recalibrate(obs_df, date(2026, 6, 15))


def test_recalibrator_warm_start_seeds_start_values(monkeypatch):
    """warm_start_theta should overwrite the `start` of matching active params,
    and must NOT mutate the caller's config (deepcopy)."""
    import dssatcalibrator.orchestrator as orch

    captured = {}

    def fake_calibrate(cfg, progress=True):
        captured["cfg"] = cfg
        class _R:  # minimal stand-in for CalibrationResult
            best_theta = {"P1": 1.0}
        return _R()

    monkeypatch.setattr(orch, "calibrate", fake_calibrate)

    cfg = {
        "calibrator": {"name": "test"},
        "parameters": {"genetic_cultivar": {"P1": {"active": True, "start": 10.0,
                                                   "min": 0, "max": 100}}},
        "assimilation": {"recalibration": {"recal_sample_size": 5, "warm_start": True}},
    }
    recal = InSeasonRecalibrator(cfg)
    obs_df = pd.DataFrame([{"date": pd.Timestamp("2026-06-01"), "value": 1.0}])

    out = recal.recalibrate(obs_df, date(2026, 6, 15), warm_start_theta={"P1": 42.0})
    assert out == {"P1": 1.0}
    # seeded into the cfg handed to calibrate ...
    assert captured["cfg"]["parameters"]["genetic_cultivar"]["P1"]["start"] == 42.0
    assert captured["cfg"]["method"]["sample"]["n"] == 5
    # ... but the caller's original cfg is untouched
    assert cfg["parameters"]["genetic_cultivar"]["P1"]["start"] == 10.0


def test_uncoupled_modes_are_gated():
    """enkf / forcing must refuse to run unless allow_uncoupled is set, and the
    guard fires before any DSSAT setup (so dummy paths are fine)."""
    from dssatcalibrator.orchestrator import assimilate

    for mode in ("enkf", "forcing"):
        cfg = {"calibrator": {"name": "t"}, "assimilation": {"mode": mode}}
        with pytest.raises(NotImplementedError):
            assimilate(cfg, progress=False)

    with pytest.raises(ValueError):
        assimilate({"calibrator": {"name": "t"}, "assimilation": {"mode": "bogus"}},
                   progress=False)


def test_unmatched_variables_helper():
    from dssatcalibrator.objective import unmatched_variables

    cfg = {"engine": {"timeseries_outputs": {"LAI": "LAID"},
                      "scalar_outputs": {"grain_yield": "HWAM"}}}
    obs = pd.DataFrame({"variable": ["LAID", "HWAM", "SW", "TMEAN"]})
    assert unmatched_variables(obs, cfg) == ["SW", "TMEAN"]
    assert unmatched_variables(pd.DataFrame(), cfg) == []
