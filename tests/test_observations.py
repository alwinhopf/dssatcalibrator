"""Tests for the observations reader against the real hemp FileA/FileT."""
import pandas as pd

from dssatcalibrator.observations import Observations, read_filea, read_filet, SCHEMA


def test_read_filea_yuku2101(hemp_dir):
    fa = read_filea(hemp_dir / "YUKU2101.HMA", "YUKU2101")
    assert list(fa.columns) == SCHEMA
    # ADAT present as a phenology event; treatment 1 anthesis = DOY 222, 2021
    adat = fa[(fa.variable == "ADAT") & (fa.treatment == 1)].iloc[0]
    assert adat["kind"] == "phenology"
    assert adat["date"] == pd.Timestamp("2021-08-10")
    # SWAH present as a scalar
    assert ((fa.variable == "SWAH") & (fa.kind == "scalar")).any()


def test_read_filea_yuba2101_doy_dates(hemp_dir):
    fa = read_filea(hemp_dir / "YUBA2101.HMA", "YUBA2101")
    edate = fa[(fa.variable == "EDATE") & (fa.treatment == 1)].iloc[0]
    adat = fa[(fa.variable == "ADAT") & (fa.treatment == 1)].iloc[0]
    assert edate["kind"] == "phenology"
    assert edate["date"] == pd.Timestamp("2021-05-17")
    assert adat["kind"] == "phenology"
    assert adat["date"] == pd.Timestamp("2021-08-04")


def test_read_filea_skips_placeholder_duplicate_headers(hemp_dir):
    fa = read_filea(hemp_dir / "UKAB2101.HMA", "UKAB2101")
    assert not fa.empty
    assert not fa["variable"].astype(str).str.startswith("__skip").any()
    assert {"ADAT", "EDAT"} <= set(fa["variable"])


def test_read_filet_yuku2101(hemp_dir):
    ft = read_filet(hemp_dir / "YUKU2101.HMT", "YUKU2101")
    assert not ft.empty
    assert (ft.kind == "timeseries").all()
    # CWAD time-series exists for several treatments with replicate rows
    cwad = ft[ft.variable == "CWAD"]
    assert set(cwad.treatment.unique()) >= {1, 2, 3}
    # replicates: more than one row per (treatment, date) for trt 1
    g = cwad[cwad.treatment == 1].groupby("date").size()
    assert (g > 1).any()


def test_observations_coverage(hemp_dir, target_experiments):
    obs = Observations.from_dssat(hemp_dir, target_experiments, crop_ext="HM")
    present = set(obs.experiments())
    # YUBA2201 ships neither FileA nor FileT -> no observations
    assert "YUBA2201" not in present
    # these have at least FileA
    for exp in ("YUFE2101", "YUFE2201", "YUKU2101", "YUKU2201", "YUBA2101"):
        assert exp in present, exp
    # CNKU2101 has only FileT
    assert "CNKU2101" in present
    cov = obs.coverage()
    assert not cov.empty
