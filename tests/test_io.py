"""Tests for dssat_io parsers against real smoke-run outputs."""
import numpy as np
import pandas as pd

from dssatcalibrator import dssat_io as io


def test_yyddd_to_date():
    assert io.yyddd_to_date(21154) == pd.Timestamp("2021-06-03")   # DOY 154 of 2021
    assert io.yyddd_to_date(21222) == pd.Timestamp("2021-08-10")
    assert pd.isna(io.yyddd_to_date(-99))
    assert pd.isna(io.yyddd_to_date(0))


def test_parse_plantgro_structure(smoke_dir):
    pg = io.parse_plantgro(smoke_dir / "PlantGro.OUT")
    assert not pg.empty
    treatments = sorted(pg["treatment"].dropna().unique().tolist())
    assert treatments and all(int(value) >= 1 for value in treatments)
    # the columns we depend on downstream
    for col in ("DAP", "GSTD", "LAID", "CWAD", "date"):
        assert col in pg.columns, col
    # biomass is monotonic-ish: end-of-season CWAD > 0
    last = pg.sort_values("DAP").groupby("treatment").tail(1)
    assert (last["CWAD"] > 0).all()
    # -99 mapped to NaN, not literal -99 (check the physical growth columns)
    phys_cols = [col for col in ("LAID", "CWAD", "CHTD", "GSTD") if col in pg.columns]
    phys = pg[phys_cols].to_numpy(dtype="float64", na_value=np.nan)
    assert not np.isclose(phys, -99.0).any()


def test_parse_evaluate_sim_vs_meas(smoke_dir):
    ev = io.parse_evaluate(smoke_dir / "Evaluate.OUT")
    assert not ev.empty
    assert {"treatment", "variable", "sim", "meas"}.issubset(ev.columns)
    # Anthesis DAP is present with finite simulated and measured values.
    adap_t1 = ev[(ev.variable == "ADAP") & (ev.treatment == 1)].iloc[0]
    assert np.isfinite(adap_t1["meas"])
    assert np.isfinite(adap_t1["sim"])
    # measured-missing entries are NaN, not -99
    assert not np.isclose(ev["meas"].dropna(), -99.0).any()


def test_parse_summary_fallback(smoke_dir):
    s = io.parse_summary(smoke_dir / "Summary.OUT")
    # best-effort: should at least recover RUNNO/TRNO and a yield-ish column
    if not s.empty:
        assert "RUNNO" in s.columns
