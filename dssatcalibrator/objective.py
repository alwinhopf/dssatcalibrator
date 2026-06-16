"""Multi-variable, multi-experiment objective: align sim vs observed and score.

Two alignment paths (see CONCEPT.md §8):

* **scalars / phenology** come straight from ``Evaluate.OUT`` (DSSAT pairs its
  own simulated & measured columns, e.g. ``ADAPS``/``ADAPM``), so no date math.
* **time-series** are matched from ``PlantGro.OUT`` to the FileT/CSV observations
  by (treatment, date), averaging replicate observations.

Scoring exposes both a minimisation ``score`` (per the configured weighting mode)
and a maximisation ``loglik`` (sigma-weighted Gaussian) so optimisers/GLUE and
the Bayesian engines share one residual table. Per-variable goodness-of-fit
metrics reuse the gridded tutorial's set: RMSE, nRMSE%, MBE, Willmott d, EF, R².
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd


def variable_maps(cfg: dict):
    ts = (cfg.get("engine", {}) or {}).get("timeseries_outputs", {}) or {}
    sc = (cfg.get("engine", {}) or {}).get("scalar_outputs", {}) or {}
    return ts, sc, {v: k for k, v in ts.items()}, {v: k for k, v in sc.items()}


def metrics(obs, sim) -> dict:
    """RMSE / nRMSE% / MBE / Willmott d / modelling efficiency EF / R² + n."""
    obs = np.asarray(obs, float)
    sim = np.asarray(sim, float)
    ok = np.isfinite(obs) & np.isfinite(sim)
    o, s = obs[ok], sim[ok]
    n = len(o)
    base = {"n": n, "RMSE": np.nan, "nRMSE_pct": np.nan, "MBE": np.nan,
            "d": np.nan, "EF": np.nan, "R2": np.nan}
    if n == 0:
        return base
    rmse = float(np.sqrt(np.mean((s - o) ** 2)))
    ob = float(np.mean(o))
    base["RMSE"] = rmse
    base["nRMSE_pct"] = 100 * rmse / ob if ob != 0 else np.nan
    base["MBE"] = float(np.mean(s - o))
    d_den = float(np.sum((np.abs(s - ob) + np.abs(o - ob)) ** 2))
    base["d"] = 1 - float(np.sum((s - o) ** 2)) / d_den if d_den > 0 else np.nan
    o_var = float(np.sum((o - ob) ** 2))
    base["EF"] = 1 - float(np.sum((o - s) ** 2)) / o_var if o_var > 0 else np.nan
    if n > 1 and np.std(o) > 0 and np.std(s) > 0:
        base["R2"] = float(np.corrcoef(o, s)[0, 1] ** 2)
    return base


def _sigma(user_var: str, obs_value: float, cfg: dict) -> float:
    spec = (cfg.get("objective", {}).get("error_model", {}) or {}).get(user_var)
    if spec is None:
        return max(abs(0.10 * obs_value), 1e-6)          # default: 10% relative
    if str(spec.get("type", "relative")) == "relative":
        return max(abs(float(spec["value"]) * obs_value), 1e-6)
    return float(spec["value"])


def build_residuals(results: dict, obs_table: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Assemble the residual table (one row per matched observation)."""
    _, _, ts_inv, sc_inv = variable_maps(cfg)
    rows = []

    # scalars / phenology from Evaluate.OUT
    for exp, res in results.items():
        ev = getattr(res, "evaluate", None)
        if ev is None or ev.empty:
            continue
        for _, r in ev.iterrows():
            base = r["variable"]
            if base not in sc_inv:
                continue
            sim, meas = r["sim"], r["meas"]
            if pd.isna(sim) or pd.isna(meas):
                if not pd.isna(meas) and pd.isna(sim):
                    # Penalise missing simulation for an observed variable by introducing a large residual
                    penalty_sim = meas + 1000.0
                    kind = "phenology" if base in ("ADAP", "EDAP", "MDAP") else "scalar"
                    rows.append(dict(exp_id=exp, treatment=int(r["treatment"]), user_var=sc_inv[base],
                                     dssat=base, kind=kind, date=pd.NaT, obs=float(meas), sim=float(penalty_sim)))
                continue
            kind = "phenology" if base in ("ADAP", "EDAP", "MDAP") else "scalar"
            rows.append(dict(exp_id=exp, treatment=int(r["treatment"]), user_var=sc_inv[base],
                             dssat=base, kind=kind, date=pd.NaT, obs=float(meas), sim=float(sim)))

    # time-series from PlantGro matched to FileT/CSV obs
    if obs_table is not None and not obs_table.empty:
        ts_obs = obs_table[obs_table["kind"] == "timeseries"]
        for exp, res in results.items():
            pg = getattr(res, "plantgro", None)
            if pg is None or pg.empty:
                continue
            o = ts_obs[ts_obs["exp_id"] == exp]
            if o.empty:
                continue
            oavg = o.groupby(["treatment", "date", "variable"], as_index=False)["value"].mean()
            for _, r in oavg.iterrows():
                col = r["variable"]
                if col not in pg.columns:
                    continue
                sub = pg[(pg["treatment"] == r["treatment"]) & (pg["date"] == r["date"])]
                if sub.empty:
                    continue
                simv = sub.iloc[0][col]
                if pd.isna(simv):
                    continue
                rows.append(dict(exp_id=exp, treatment=int(r["treatment"]),
                                 user_var=ts_inv.get(col, col), dssat=col, kind="timeseries",
                                 date=r["date"], obs=float(r["value"]), sim=float(simv)))

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["sigma"] = [_sigma(uv, ov, cfg) for uv, ov in zip(df["user_var"], df["obs"])]
    wts = cfg.get("objective", {}).get("weights", {}) or {}
    df["weight"] = df["user_var"].map(lambda v: float(wts.get(v, 1.0)))
    df["resid"] = df["sim"] - df["obs"]
    if cfg.get("objective", {}).get("obs_autocorr", False):
        df = _downweight_autocorr(df)
    return df


def _downweight_autocorr(df: pd.DataFrame) -> pd.DataFrame:
    """Down-weight dense time-series for serial correlation (``obs_autocorr: true``).

    Consecutive daily LAI/biomass points are *not* independent measurements: if the
    model is too high on day 50 it is almost certainly too high on day 51. Treating
    every day as a fresh observation makes a long time-series dominate one-off
    scalars like grain yield. For each ``(experiment, variable, treatment)`` series
    we estimate a lag-1 autocorrelation ``rho`` and shrink every point's weight by
    the AR(1) effective-sample-size factor ``(1 - rho) / (1 + rho)`` — so a strongly
    correlated 100-point series counts for only a handful of *effective* points.
    """
    df = df.copy()
    ts = df[df["kind"] == "timeseries"]
    for (_exp, _uv, _trt), g in ts.groupby(["exp_id", "user_var", "treatment"]):
        if len(g) < 3:
            continue
        x = g.sort_values("date")["obs"].to_numpy(dtype=float)
        x = x - x.mean()
        denom = float(np.sum(x * x))
        if denom <= 0:
            continue
        rho = float(np.sum(x[:-1] * x[1:]) / denom)
        rho = min(max(rho, 0.0), 0.99)          # only down-weight positive correlation
        factor = (1.0 - rho) / (1.0 + rho)
        df.loc[g.index, "weight"] = df.loc[g.index, "weight"] * factor
    return df


@dataclass
class ObjectiveResult:
    score: float                       # to MINIMISE
    loglik: float                      # to MAXIMISE (sigma-weighted Gaussian)
    residuals: pd.DataFrame
    per_var: dict = field(default_factory=dict)
    per_exp_var: pd.DataFrame = field(default_factory=pd.DataFrame)


def score(results: dict, obs_table: pd.DataFrame, cfg: dict) -> ObjectiveResult:
    """Score a set of per-experiment spawn results against observations."""
    resid = build_residuals(results, obs_table, cfg)
    if resid.empty:
        return ObjectiveResult(score=float("inf"), loglik=float("-inf"),
                               residuals=resid, per_var={}, per_exp_var=pd.DataFrame())

    weighting = cfg.get("objective", {}).get("weighting", "unified")
    wts = cfg.get("objective", {}).get("weights", {}) or {}

    loglik = float(-0.5 * np.sum(((resid["resid"] / resid["sigma"]) ** 2) * resid["weight"]))

    if weighting == "sigma":
        # Total sigma-weighted chi-square (NO per-variable averaging). The
        # statistically 'correct' misfit; matches the Gaussian log-likelihood, and
        # is what the Bayesian engines expect.
        sc = float(np.sum(((resid["resid"] / resid["sigma"]) ** 2) * resid["weight"]))
    elif weighting == "user":
        # Raw normalised RMSE per variable x explicit weights. Full manual control,
        # no automatic balancing of counts or scales.
        sc = 0.0
        for uv, g in resid.groupby("user_var"):
            ob = g["obs"].mean()
            nrmse = float(np.sqrt(np.mean(g["resid"] ** 2))) / (abs(ob) if ob else 1.0)
            sc += float(wts.get(uv, 1.0)) * nrmse
    elif weighting == "count_scale":
        # Each variable contributes the AVERAGE normalised squared error of its own
        # points (count-balanced: a 100-point series can't drown a 1-point yield),
        # then we average ACROSS variables -> a clean 'mean per-variable misfit'.
        per = [float(wts.get(uv, 1.0)) * float(np.mean((g["resid"] / g["sigma"]) ** 2))
               for uv, g in resid.groupby("user_var")]
        sc = float(np.mean(per)) if per else float("inf")
    else:  # "unified" (default) and "agmip_wls"
        # SUM of per-variable mean normalised squared errors x weights. Count- and
        # scale-balanced like count_scale, but summed (so fitting more variables
        # costs more). 'agmip_wls' uses this same surface after the driver has reset
        # the per-variable weights to 1 / residual-variance (see pipeline).
        sc = 0.0
        for uv, g in resid.groupby("user_var"):
            mse = float(np.mean((g["resid"] / g["sigma"]) ** 2))
            sc += float(wts.get(uv, 1.0)) * mse

    per_var = {uv: metrics(g["obs"], g["sim"]) for uv, g in resid.groupby("user_var")}
    pev = (resid.groupby(["exp_id", "user_var"])
           .apply(lambda g: pd.Series(metrics(g["obs"], g["sim"])), include_groups=False)
           .reset_index())
    return ObjectiveResult(score=float(sc), loglik=loglik, residuals=resid,
                           per_var=per_var, per_exp_var=pev)
