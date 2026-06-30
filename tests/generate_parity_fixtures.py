"""Generate cross-language parity fixtures from the Python implementation.

The Python package is the source of truth. This writes JSON golden files that
the R testthat suite (tests/testthat/) loads and asserts its own output against,
so R↔Python parity is checked on the user's machine where R is available.

Run from the repo root:  PYTHONPATH=python python tests/generate_parity_fixtures.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from dssatcalibrator import priors, objective, config as cfgmod
from dssatcalibrator.spaces import ParameterSpace

OUT = Path(__file__).parent / "fixtures"
OUT.mkdir(exist_ok=True)


def _sanitize(o):
    """Replace non-finite floats with string sentinels so the JSON is standard
    (R's jsonlite rejects bare Infinity/NaN). The R tests decode these back."""
    if isinstance(o, (float, np.floating)):
        o = float(o)
        if o != o:
            return "nan"
        if o == float("inf"):
            return "inf"
        if o == float("-inf"):
            return "-inf"
        return float(f"{o:.15g}")
    if isinstance(o, dict):
        return {k: _sanitize(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_sanitize(v) for v in o]
    return o


def dump(name: str, obj) -> None:
    with open(OUT / name, "w", encoding="utf-8") as fh:
        json.dump(_sanitize(obj), fh, indent=2, default=str, allow_nan=False)
    print("wrote", OUT / name)


# --- priors.log_prior_one (deterministic; exact parity expected) -------------
prior_specs = {
    "uniform":     {"name": "U", "min": 0.0,  "max": 10.0},
    "normal":      {"name": "N", "min": 0.0,  "max": 10.0, "start": 5.0,
                    "prior": {"dist": "normal", "sd": 2.0}},
    "lognormal":   {"name": "L", "min": 0.1,  "max": 50.0, "start": 5.0,
                    "prior": {"dist": "lognormal", "sd": 0.5}},
    "triangular":  {"name": "T", "min": 0.0,  "max": 10.0, "start": 7.0,
                    "prior": {"dist": "triangular"}},
}
test_points = [-1.0, 0.05, 0.5, 1.0, 3.0, 5.0, 7.0, 9.5, 11.0]
prior_gold = {}
for dist, spec in prior_specs.items():
    prior_gold[dist] = {
        "spec": spec,
        "points": test_points,
        "log_prior": [priors.log_prior_one(spec, v) for v in test_points],
    }
dump("priors_logdensity.json", prior_gold)


# --- objective.metrics (pure; exact parity expected) -------------------------
cases = {
    "perfect":   {"obs": [1, 2, 3, 4, 5],        "sim": [1, 2, 3, 4, 5]},
    "biased":    {"obs": [10, 20, 30, 40],       "sim": [12, 22, 33, 39]},
    "noisy":     {"obs": [5, 7, 9, 11, 13, 15],  "sim": [4, 8, 8, 12, 12, 16]},
    "single":    {"obs": [42.0],                 "sim": [40.0]},
    "with_nan":  {"obs": [1, 2, float("nan"), 4],"sim": [1.1, 2.2, 3.0, 3.8]},
}
metrics_gold = {k: {"input": v, "metrics": objective.metrics(v["obs"], v["sim"])}
                for k, v in cases.items()}
dump("metrics.json", metrics_gold)


# --- config.load_config + active_parameters + ParameterSpace -----------------
sample_yaml = """
calibrator:
  name: parity_test
  seed: 7
method:
  preset: C
parameters:
  cultivar:
    P1:  { active: true,  min: 100, max: 500, start: 300 }
    P5:  { active: true,  min: 300, max: 700, start: 505, prior: {dist: normal, sd: 50} }
    G1:  { active: false, min: 10,  max: 40,  start: 25 }
  ecotype:
    PHINT: { active: true, min: 60, max: 120 }
"""
cfg_path = OUT / "_sample_config.yaml"
cfg_path.write_text(sample_yaml, encoding="utf-8")
cfg = cfgmod.load_config(cfg_path)
act = cfgmod.active_parameters(cfg)
space = ParameterSpace.from_config(cfg)
dump("config_space.json", {
    "active_names": [a["name"] for a in act],
    "active_groups": [a["group"] for a in act],
    "space_names": list(space.names),
    "space_low": [float(x) for x in space.low],
    "space_high": [float(x) for x in space.high],
    "space_start": [float(x) for x in space.start],
    "ndim": space.ndim,
    "merged_preset": cfg["method"]["preset"],
    "merged_seed": cfg["calibrator"]["seed"],
    # a default that must survive the merge untouched:
    "default_gating_ecotype": cfg["gating"]["ecotype"],
})

# --- Phase 2: dssat_io.yyddd_to_date (deterministic) -------------------------
from dssatcalibrator.dssat_io import yyddd_to_date
date_codes = [-99, 0, 98032, 1, 24001, 1998032, 2024100, 99366, 50180, 80001]
def _isodate(ts):
    return None if (ts is None or (hasattr(ts, "__class__") and str(ts) == "NaT")) else ts.strftime("%Y-%m-%d")
import pandas as _pd
dates_gold = {}
for c in date_codes:
    d = yyddd_to_date(c)
    dates_gold[str(c)] = None if _pd.isna(d) else d.strftime("%Y-%m-%d")
dump("yyddd_dates.json", dates_gold)


# --- Phase 2: per-source error_model + obs operator (deterministic) ----------
from dssatcalibrator.sources.field import FieldMeasurementSource
from dssatcalibrator.sources.uav import UAVMultispectralSource
from dssatcalibrator.sources.iot import SoilMoistureSensorSource, CanopyTemperatureSource
from dssatcalibrator.sources.satellite import (SentinelLAISource, MODISLAISource,
                                               _apply_obs_operator)
from dssatcalibrator.sources.farm_software import FarmPhenologySource, FarmManagementSource

err_gold = {}
def errcase(label, src, var, val, meta):
    err_gold[label] = {"variable": var, "value": val, "metadata": meta,
                       "sigma": src.error_model(var, val, meta)}

f = FieldMeasurementSource({})
errcase("field_LAID", f, "LAID", 3.2, {})
errcase("field_HWAM", f, "HWAM", 4000.0, {})
errcase("field_ADAT", f, "ADAT", 99001.0, {})
errcase("field_unknown", f, "ZZZZ", 50.0, {})

u = UAVMultispectralSource({})
errcase("uav_LAID_good", u, "LAID", 2.0, {"flight_quality": "good"})
errcase("uav_LAID_poor", u, "LAID", 2.0, {"flight_quality": "poor"})
errcase("uav_unknown", u, "WWWW", 10.0, {"flight_quality": "good"})

sm = SoilMoistureSensorSource({})
errcase("iot_sw_cap_factory", sm, "SW", 0.3, {"sensor_type": "capacitance", "calibration_status": "factory"})
errcase("iot_sw_tdr_calib", sm, "SW", 0.3, {"sensor_type": "tdr", "calibration_status": "field_calibrated"})

ct = CanopyTemperatureSource({})
errcase("iot_tmean", ct, "TMEAN", 25.0, {})

s2 = SentinelLAISource({})
errcase("sentinel_below_sat", s2, "LAID", 2.0, {"cloud_fraction": 0.0})
errcase("sentinel_above_sat_cloud", s2, "LAID", 6.0, {"cloud_fraction": 0.4})

mo = MODISLAISource({})
errcase("modis_qc0", mo, "LAID", 3.0, {"qc_flag": 0})
errcase("modis_qc2", mo, "LAID", 3.0, {"qc_flag": 2})

fp = FarmPhenologySource({})
errcase("farm_gstd_weekly", fp, "GSTD", 5.0, {"date_precision": "weekly"})
errcase("farm_adat_biweekly", fp, "ADAT", 1.0, {"date_precision": "biweekly"})

fm = FarmManagementSource({})
errcase("mgmt_date", fm, "spray_date", 1.0, {})
errcase("mgmt_amount", fm, "irrig_amount", 40.0, {})

err_gold["_obs_operator"] = {
    "identity": _apply_obs_operator({}, 3.0),
    "scaled":   _apply_obs_operator({"obs_operator": {"scale": 0.9, "offset": 0.2}}, 3.0),
}
dump("source_error_models.json", err_gold)


# --- Phase 2: inverse-variance fusion merge (deterministic) ------------------
from dssatcalibrator.fusion import ObservationFuser
import pandas as pd
fuse_df = pd.DataFrame([
    # two coincident LAID obs on same exp/trt/date -> merged by inverse variance
    dict(exp_id="E1", treatment=1, variable="LAID", kind="timeseries",
         date=pd.Timestamp("2021-06-01"), value=2.0, sigma=0.5, weight=1.0,
         source="sentinel2_lai", quality_flag=0, spatial_res_m=10.0),
    dict(exp_id="E1", treatment=1, variable="LAID", kind="timeseries",
         date=pd.Timestamp("2021-06-01"), value=2.6, sigma=0.7, weight=1.0,
         source="modis_lai", quality_flag=0, spatial_res_m=250.0),
    # a lone obs -> passes through unchanged
    dict(exp_id="E1", treatment=1, variable="CWAD", kind="timeseries",
         date=pd.Timestamp("2021-06-10"), value=1500.0, sigma=120.0, weight=1.0,
         source="field_measurements", quality_flag=0, spatial_res_m=float("nan")),
])
fuser = ObservationFuser([], {"fusion": {"conflict_resolution": "inverse_variance"}})
merged = fuser._inverse_variance_merge(fuse_df)
merged = merged.sort_values(["variable", "date"]).reset_index(drop=True)
dump("fusion_inverse_variance.json", {
    "input": [{**{k: (v.strftime("%Y-%m-%d") if isinstance(v, pd.Timestamp) else v)
                  for k, v in row.items()}} for row in fuse_df.to_dict("records")],
    "merged": [{"variable": r["variable"], "date": r["date"].strftime("%Y-%m-%d"),
                "value": float(r["value"]), "sigma": float(r["sigma"]),
                "source": r["source"]} for _, r in merged.iterrows()],
})


# --- Phase 3: writers (byte-level fixed-width editing) -----------------------
from dssatcalibrator import writers
import shutil, tempfile

# Build column-aligned DSSAT rows from the header (values end-justified to each
# header token's end column, anchor code left-justified at column 0) — this is
# how DSSAT writes fixed-width files, so the editors land cleanly.
def _aligned_row(header, anchor, coeffs, eco_code=None, vrname=None):
    bounds = writers.parse_header_boundaries(header)
    width = max(hi for _, hi in bounds.values())
    chars = [" "] * width
    chars[0:len(anchor)] = list(anchor)                       # VAR#/ECO# code at col 0
    if vrname:                                                # name just after code + space
        nlo = len(anchor) + 1
        chars[nlo:nlo + len(vrname)] = list(vrname)
    if eco_code and "ECO#" in bounds:
        lo, hi = bounds["ECO#"]
        chars[lo:hi] = list(eco_code.rjust(hi - lo)[:hi - lo])
    for name, val in coeffs.items():
        lo, hi = bounds[name]
        chars[lo:hi] = list(writers._fmt(float(val), hi - lo))
    return "".join(chars).rstrip()

CUL_HDR = "@VAR#  VRNAME.......... EXPNO   ECO#  CSDL PPSEN EM-FL FL-SH FL-SD SD-PM FL-LF LFMAX"
_cul_coeffs = {"CSDL": 12.33, "PPSEN": 0.249, "EM-FL": 18.50, "FL-SH": 6.0,
               "FL-SD": 14.0, "SD-PM": 32.0, "FL-LF": 18.0, "LFMAX": 1.0}
_cul_min = {"CSDL": 11.0, "PPSEN": 0.1, "EM-FL": 15.0, "FL-SH": 4.0,
            "FL-SD": 10.0, "SD-PM": 28.0, "FL-LF": 14.0, "LFMAX": 0.8}
_cul_max = {"CSDL": 14.0, "PPSEN": 0.4, "EM-FL": 22.0, "FL-SH": 8.0,
            "FL-SD": 18.0, "SD-PM": 36.0, "FL-LF": 22.0, "LFMAX": 1.4}
SAMPLE_CUL = "\n".join([
    "*CULTIVARS:CRGRO048", "! sample for parity tests", CUL_HDR,
    _aligned_row(CUL_HDR, "IB0008", _cul_coeffs, "SB0301", "YUNMA8"),
    _aligned_row(CUL_HDR, "999991", _cul_min, "SB0301", "MINIMA"),
    _aligned_row(CUL_HDR, "999992", _cul_max, "SB0301", "MAXIMA"),
]) + "\n"

ECO_HDR = "@ECO#  ECONAME.......... PL-EM EM-V1 V1-JU JU-R0  PM06  PM09"
_eco_coeffs = {"PL-EM": 3.6, "EM-V1": 3.0, "V1-JU": 6.0, "JU-R0": 12.0,
               "PM06": 0.0, "PM09": 0.0}
SAMPLE_ECO = "\n".join([
    "*ECOTYPE:CRGRO048", ECO_HDR,
    _aligned_row(ECO_HDR, "SB0301", _eco_coeffs, None, "GENERIC SOYBEAN"),
]) + "\n"

(OUT / "sample.CUL").write_text(SAMPLE_CUL, encoding="utf-8")
(OUT / "sample.ECO").write_text(SAMPLE_ECO, encoding="utf-8")

def _edit_after(sample_name, edit_fn, anchor, updates):
    tmpd = Path(tempfile.mkdtemp())
    tmp = tmpd / sample_name
    shutil.copy(OUT / sample_name, tmp)
    edit_fn(tmp, anchor, updates)
    return tmp.read_text(encoding="utf-8")

cul_updates = {"CSDL": 12.5, "LFMAX": 1.05, "EM-FL": 19.0}
cul_after = _edit_after("sample.CUL", writers.edit_cultivar, "IB0008", cul_updates)
(OUT / "sample_cul_after.CUL").write_text(cul_after, encoding="utf-8")

eco_updates = {"PL-EM": 3.8, "JU-R0": 11.5}
eco_after = _edit_after("sample.ECO", writers.edit_ecotype, "SB0301", eco_updates)
(OUT / "sample_eco_after.ECO").write_text(eco_after, encoding="utf-8")

# field maps, bounds, read-back, _fmt, header boundaries
fmap_cul = writers.cultivar_field_map(OUT / "sample.CUL")
bounds = writers.read_cul_calibration_bounds(OUT / "sample.CUL")
read_back = writers.read_cultivar_values(OUT / "sample.CUL", "IB0008")
hdr = "@P    PDATE EDATE  PPOP  PPOE  PLME  PLDS  PLRS  PLRD  PLDP"
hbound = writers.parse_header_boundaries(hdr)
fmt_cases = {f"{v}_{w}": writers._fmt(v, w)
             for v, w in [(12.5, 5), (1.05, 5), (1234.5, 4), (0.249, 5),
                          (-99.0, 6), (3.0, 6), (18.0, 5), (0.18, 5)]}
dump("writers_meta.json", {
    "cul_updates": cul_updates,
    "eco_updates": eco_updates,
    "cultivar_field_map": {k: list(v) for k, v in fmap_cul.items()},
    "calibration_bounds": bounds,
    "read_cultivar_values": {k: v for k, v in read_back.items() if v is not None},
    "parse_header_boundaries": {k: list(v) for k, v in hbound.items()},
    "fmt_cases": fmt_cases,
})


# --- Phase 3b: dssat_io parsers + spawn helpers ------------------------------
from dssatcalibrator import dssat_io, spawn

SAMPLE_PLANTGRO = (
"*RUN   1        : parity sample\n"
" MODEL          CRGRO048\n"
" TREATMENT  1   RAINFED\n"
"\n"
"@YEAR DOY   DAS   DAP  L#SD  GSTD  LAID  CWAD  GWAD\n"
" 2021 152     0     0   0.0   1.0  0.00     0     0\n"
" 2021 160     8     8   2.0   2.0  0.45   120     0\n"
" 2021 200    48    48   8.0   5.0  3.20  3500   100\n"
)
SAMPLE_EVALUATE = (
"@RUN EXCODE        TN RN CR ADAPS ADAPM HWAMS HWAMM LAIXS LAIXM\n"
"   1 SAMPLE0001     1  1 SB    45    47  3200  3100  3.20  3.00\n"
"   2 SAMPLE0001     2  1 SB    50    49  3400  3500  3.50  3.40\n"
)
SAMPLE_SUMMARY = (
"@  RUNNO   TRNO    CR  TNAM    CWAM   HWAM   ADAT   MDAT\n"
"      1      1    SB  TRT1    3500   3200  21045  21120\n"
"      2      2    SB  TRT2    3600   3400  21050  21125\n"
)
(OUT / "PlantGro.OUT").write_text(SAMPLE_PLANTGRO, encoding="utf-8")
(OUT / "Evaluate.OUT").write_text(SAMPLE_EVALUATE, encoding="utf-8")
(OUT / "Summary.OUT").write_text(SAMPLE_SUMMARY, encoding="utf-8")

pg = dssat_io.parse_plantgro(OUT / "PlantGro.OUT")
ev = dssat_io.parse_evaluate(OUT / "Evaluate.OUT")
sm = dssat_io.parse_summary(OUT / "Summary.OUT")

dump("parsers.json", {
    "plantgro": {
        "nrow": int(len(pg)),
        "columns": list(pg.columns),
        "treatment": [int(x) for x in pg["treatment"]],
        "LAID": [float(x) for x in pg["LAID"]],
        "CWAD": [float(x) for x in pg["CWAD"]],
        "dates": [d.strftime("%Y-%m-%d") for d in pg["date"]],
    },
    "evaluate": [
        {"treatment": int(r["treatment"]), "variable": r["variable"],
         "sim": float(r["sim"]), "meas": float(r["meas"])}
        for _, r in ev.sort_values(["treatment", "variable"]).iterrows()
    ],
    "summary": {
        "nrow": int(len(sm)),
        "RUNNO": [int(x) for x in sm["RUNNO"]],
        "CWAM": [float(x) for x in sm["CWAM"]],
        "HWAM": [float(x) for x in sm["HWAM"]],
    },
})

# spawn helpers
thetas = {
    "t1": {"P1": 380, "P5": 505.0, "G1": 26.5},
    "t2": {"CSDL": 12.33, "LFMAX": 1.0, "EM-FL": 18.5},
    "t3": {"x": 0.249, "y": 1000.0},
}
# write_dssbatch output (read the file it writes)
tmpd = Path(tempfile.mkdtemp())
batch = spawn.write_dssbatch(tmpd, "EXP0001.HMX", [1, 2, 10])
dump("spawn_helpers.json", {
    "theta_hash": {k: spawn.theta_hash(v) for k, v in thetas.items()},
    "thetas": thetas,
    "write_dssbatch": batch.read_text(),
    "normalize": {
        "dedup_order": spawn._normalize_treatments([3, 1, 1, 10, 3], "native"),
        "single": spawn._normalize_treatments([5], "native"),
    },
})


# --- Phase 4: GLUE post-processing (deterministic) ---------------------------
from dssatcalibrator.engines import run_glue
import numpy as _np
design_glue = pd.DataFrame({
    "sample_id": [0, 1, 2, 3, 4],
    "P1": [300.0, 350.0, 380.0, 410.0, 500.0],
    "P5": [450.0, 500.0, 505.0, 520.0, 600.0],
    "score": [2.5, 1.2, 0.8, 1.0, _np.inf],
    "loglik": [-12.5, -6.0, -4.0, -5.0, -_np.inf],
})
g = run_glue(design_glue.copy(), ["P1", "P5"], {"method": {"bayesian": {"behavioural_quantile": 0.4}}}, space=None)
dump("glue.json", {
    "input": {k: [None if (isinstance(v, float) and not _np.isfinite(v)) else v for v in design_glue[k].tolist()]
              for k in design_glue.columns},
    "behavioural_quantile": 0.4,
    "weight": [float(x) for x in g.design["weight"]],
    "ess": float(g.ess),
    "threshold": float(g.threshold),
    "best_sample_id": int(g.best_sample_id),
    "best_theta": {k: float(v) for k, v in g.best_theta.items()},
    "n_behavioural": int(len(g.behavioural)),
})


print("\nAll parity fixtures generated.")
