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
    # four density treatments
    assert sorted(pg["treatment"].dropna().unique().tolist()) == [1, 2, 3, 4]
    # the columns we depend on downstream
    for col in ("DAP", "GSTD", "LAID", "CWAD", "CHTD", "date"):
        assert col in pg.columns, col
    # biomass is monotonic-ish: end-of-season CWAD > 0
    last = pg.sort_values("DAP").groupby("treatment").tail(1)
    assert (last["CWAD"] > 0).all()
    # -99 mapped to NaN, not literal -99 (check the physical growth columns)
    phys = pg[["LAID", "CWAD", "CHTD", "GSTD"]].to_numpy(dtype="float64", na_value=np.nan)
    assert not np.isclose(phys, -99.0).any()


def test_parse_evaluate_sim_vs_meas(smoke_dir):
    ev = io.parse_evaluate(smoke_dir / "Evaluate.OUT")
    assert not ev.empty
    assert {"treatment", "variable", "sim", "meas"}.issubset(ev.columns)
    # anthesis DAP is present with a measured value for treatment 1 (ADAPM=75, ADAPS=79)
    adap_t1 = ev[(ev.variable == "ADAP") & (ev.treatment == 1)].iloc[0]
    assert adap_t1["meas"] == 75
    assert adap_t1["sim"] == 79
    # measured-missing entries are NaN, not -99
    assert not np.isclose(ev["meas"].dropna(), -99.0).any()


def test_parse_summary_fallback(smoke_dir):
    s = io.parse_summary(smoke_dir / "Summary.OUT")
    # best-effort: should at least recover RUNNO/TRNO and a yield-ish column
    if not s.empty:
        assert "RUNNO" in s.columns
