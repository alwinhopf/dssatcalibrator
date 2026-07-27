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


_PHENOLOGY_DATE_TO_DAP = {
    "EDAT": "EDAP",
    "EDATE": "EDAP",
    "ADAT": "ADAP",
    "MDAT": "MDAP",
    "HDAT": "HDAP",
    "R8": "R8AP",
}
_PHENOLOGY_DAP_OUTPUTS = {"EDAP", "ADAP", "MDAP", "HDAP", "R8AP"}


def unmatched_variables(obs_table: pd.DataFrame, cfg: dict) -> list[str]:
    """DSSAT variables present in the observations but NOT scorable.

    A variable can only be scored if it appears in ``engine.timeseries_outputs``
    or ``engine.scalar_outputs`` (and is produced by the spawn parser). Sources
    such as the IoT/UAV adapters can ingest variables (``SW``, ``TMEAN``,
    ``canopy_cover`` …) that no configured output maps to; those rows are silently
    dropped when scoring. This returns them so the orchestrator can warn the user.
    """
    ts, sc, _, _ = variable_maps(cfg)
    known = set(ts.values()) | set(sc.values())
    if obs_table is None or getattr(obs_table, "empty", True) or "variable" not in obs_table:
        return []
    present = {str(v) for v in obs_table["variable"].dropna().unique()}
    mapped = {
        v for v in present
        if _PHENOLOGY_DATE_TO_DAP.get(v) in known
    }
    return sorted(present - known - mapped)


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
        sigma = max(abs(0.10 * obs_value), 1e-6)          # default: 10% relative
    elif str(spec.get("type", "relative")) == "relative":
        sigma = max(abs(float(spec["value"]) * obs_value), 1e-6)
    else:
        sigma = float(spec["value"])

    # Optional model discrepancy: observations can be precise while the selected
    # DSSAT module is still structurally imperfect, especially for new species.
    # Add this in quadrature so calibration does not over-bend parameters to
    # explain known model error.
    disc_cfg = (cfg.get("objective", {}) or {}).get("model_discrepancy", {}) or {}
    disc = disc_cfg.get("default", disc_cfg.get("value", 0.0))
    variables = disc_cfg.get("variables", {}) or {}
    relative = disc_cfg.get("relative", {}) or {}
    if user_var in variables:
        disc = variables[user_var]
    if user_var in relative:
        disc = max(float(disc), abs(float(relative[user_var]) * obs_value))
    disc = float(disc or 0.0)
    if disc > 0:
        sigma = float(np.sqrt(sigma ** 2 + disc ** 2))
    return max(float(sigma), 1e-12)


def _standardized_loss(z, cfg: dict):
    """Robust squared-error replacement in standardized residual units.

    objective.likelihood:
      gaussian  -> z^2
      student_t -> heavy-tailed pseudo deviance
      huber     -> quadratic near zero, linear in the tails
    """
    lcfg = (cfg.get("objective", {}) or {}).get("likelihood", {}) or {}
    if isinstance(lcfg, str):
        kind, lcfg = lcfg.lower(), {}
    else:
        kind = str(lcfg.get("type", "gaussian")).lower()
    z = np.asarray(z, dtype=float)
    if kind in ("student_t", "student-t", "t"):
        nu = max(float(lcfg.get("df", lcfg.get("nu", 4.0))), 1.01)
        return (nu + 1.0) * np.log1p((z ** 2) / nu)
    if kind == "huber":
        delta = max(float(lcfg.get("delta", 1.5)), 1e-9)
        az = np.abs(z)
        return np.where(az <= delta, z ** 2, 2.0 * delta * az - delta ** 2)
    return z ** 2


def _weighted_loss(resid: pd.DataFrame, cfg: dict) -> pd.Series:
    z = resid["resid"].to_numpy(dtype=float) / resid["sigma"].to_numpy(dtype=float)
    return pd.Series(_standardized_loss(z, cfg), index=resid.index)


def _group_loss(group: pd.DataFrame, cfg: dict) -> float:
    """Aggregate standardized point losses for one output variable.

    ``score_metric: rmse`` is useful for calibration objectives stated directly
    as RMSE.  The default remains mean squared standardized error for backwards
    compatibility with existing configurations and Bayesian likelihoods.
    """
    loss = group["_loss"].to_numpy(dtype=float)
    weights = group.get("weight", pd.Series(1.0, index=group.index)).to_numpy(dtype=float)
    valid = np.isfinite(loss) & np.isfinite(weights) & (weights > 0)
    if valid.any():
        mean_loss = float(np.average(loss[valid], weights=weights[valid]))
    else:
        mean_loss = float(np.mean(loss))
    metric = str((cfg.get("objective", {}) or {}).get("score_metric", "mse")).lower()
    if metric in ("rmse", "root_mean_square", "root_mean_squared_error"):
        return float(np.sqrt(max(mean_loss, 0.0)))
    return mean_loss


def _drop_configured_zero_observations(resid: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Drop exact-zero observations for configured variables before scoring.

    Some FileT variables use ``0`` as a practical placeholder after the measured
    quantity has become unavailable. That is especially toxic for relative error
    models because sigma collapses toward zero. Keep this opt-in by variable so
    valid zeros, such as early grain mass, remain usable.
    """
    ocfg = (cfg.get("objective", {}) or {})
    raw = ocfg.get("ignore_zero_observations", ocfg.get("drop_zero_observations", []))
    if not raw:
        return resid
    drop_all = False
    atol = 1e-12
    if isinstance(raw, dict):
        atol = float(raw.get("atol", raw.get("tolerance", atol)))
        raw = raw.get("variables", raw.get("user_vars", raw.get("dssat_vars", [])))
    elif isinstance(raw, bool):
        drop_all = bool(raw)
        raw = []
    names = {str(v) for v in raw}
    obs = pd.to_numeric(resid["obs"], errors="coerce")
    zero = np.isclose(obs, 0.0, atol=atol, rtol=0.0)
    if drop_all:
        keep = ~zero
    else:
        variable_match = resid["user_var"].astype(str).isin(names) | resid["dssat"].astype(str).isin(names)
        keep = ~(zero & variable_match)
    return resid.loc[keep].copy()


def _planting_date_from_plantgro(pg: pd.DataFrame, treatment: int):
    if pg is None or pg.empty or "date" not in pg.columns or "DAP" not in pg.columns:
        return pd.NaT
    sub = pg[pd.to_numeric(pg.get("treatment"), errors="coerce") == int(treatment)]
    if sub.empty:
        return pd.NaT
    dates = pd.to_datetime(sub["date"], errors="coerce")
    dap = pd.to_numeric(sub["DAP"], errors="coerce")
    ok = dates.notna() & dap.notna()
    if not ok.any():
        return pd.NaT
    planted = dates[ok] - pd.to_timedelta(dap[ok], unit="D")
    return planted.iloc[0] if not planted.empty else pd.NaT


def _scalar_observation_mapping(obs_var: str, sc_map: dict, sc_inv: dict) -> tuple[str, str]:
    """Map observed FileA variable names such as ADAT onto ADAP-style outputs."""
    obs_var = str(obs_var)
    if obs_var in sc_map:
        dssat_var = sc_map[obs_var]
        return obs_var, dssat_var
    mapped = _PHENOLOGY_DATE_TO_DAP.get(obs_var)
    if mapped and mapped in sc_inv:
        return sc_inv[mapped], mapped
    if obs_var in sc_inv:
        return sc_inv[obs_var], obs_var
    return obs_var, obs_var


def _observed_scalar_value(row: pd.Series, dssat_var: str, pg: pd.DataFrame) -> float:
    if str(row.get("kind")) == "phenology" and dssat_var in _PHENOLOGY_DAP_OUTPUTS:
        date = pd.to_datetime(row.get("date"), errors="coerce")
        planted = _planting_date_from_plantgro(pg, int(row["treatment"]))
        if pd.notna(date) and pd.notna(planted):
            return float((date.normalize() - planted.normalize()).days)
    return float(row["value"])


def _timeseries_sim_value(pg: pd.DataFrame, treatment: int, date, col: str):
    sub = pg[pd.to_numeric(pg["treatment"], errors="coerce") == int(treatment)].copy()
    if sub.empty or col not in sub.columns:
        return np.nan
    dates = pd.to_datetime(sub["date"], errors="coerce")
    target = pd.Timestamp(date)
    exact = sub[dates == target]
    if not exact.empty:
        return exact.iloc[0][col]
    ok = dates.notna()
    if ok.any() and target > dates[ok].max():
        tail = sub.loc[ok].assign(_date=dates[ok]).sort_values("_date").iloc[-1]
        return tail[col]
    return np.nan


def build_residuals(results: dict, obs_table: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Assemble the residual table (one row per matched observation)."""
    ts_map, sc_map, ts_inv, sc_inv = variable_maps(cfg)
    rows = []
    seen_scalar = set()

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
                    seen_scalar.add((exp, int(r["treatment"]), sc_inv[base], kind))
                continue
            kind = "phenology" if base in ("ADAP", "EDAP", "MDAP") else "scalar"
            rows.append(dict(exp_id=exp, treatment=int(r["treatment"]), user_var=sc_inv[base],
                             dssat=base, kind=kind, date=pd.NaT, obs=float(meas), sim=float(sim)))
            seen_scalar.add((exp, int(r["treatment"]), sc_inv[base], kind))

    # CSV / fused scalar observations: match observed scalar rows to Evaluate.OUT
    # simulated values even when Evaluate.OUT does not carry measured FileA data.
    if obs_table is not None and not obs_table.empty:
        scalar_obs = obs_table[obs_table["kind"].isin(["scalar", "phenology"])]
        for exp, res in results.items():
            ev = getattr(res, "evaluate", None)
            if ev is None:
                ev = pd.DataFrame()
            o = scalar_obs[scalar_obs["exp_id"] == exp]
            if o.empty:
                continue
            pg = getattr(res, "plantgro", None)
            computed = []
            for _, rr in o.iterrows():
                user_var, dssat_var = _scalar_observation_mapping(rr["variable"], sc_map, sc_inv)
                computed.append({
                    "treatment": rr["treatment"],
                    "variable": user_var,
                    "dssat": dssat_var,
                    "kind": rr["kind"],
                    "value": _observed_scalar_value(rr, dssat_var, pg),
                })
            oavg = pd.DataFrame(computed)
            if oavg.empty:
                continue
            oavg = oavg.dropna(subset=["value"]).groupby(
                ["treatment", "variable", "dssat", "kind"], as_index=False
            )["value"].mean()
            for _, r in oavg.iterrows():
                user_var = r["variable"]
                dssat_var = r["dssat"]
                if dssat_var not in sc_inv:
                    continue
                available = set(ev["variable"]) if "variable" in ev else set()
                if dssat_var not in available:
                    if pd.isna(r["value"]):
                        continue
                    key = (exp, int(r["treatment"]), user_var, r["kind"])
                    if key in seen_scalar:
                        continue
                    rows.append(dict(exp_id=exp, treatment=int(r["treatment"]),
                                     user_var=user_var, dssat=dssat_var, kind=r["kind"],
                                     date=pd.NaT, obs=float(r["value"]),
                                     sim=float(r["value"]) + 1000.0))
                    seen_scalar.add(key)
                    continue
                key = (exp, int(r["treatment"]), user_var, r["kind"])
                if key in seen_scalar:
                    continue
                sub = ev[(ev["treatment"] == r["treatment"]) & (ev["variable"] == dssat_var)]
                if pd.isna(r["value"]):
                    continue
                if sub.empty or pd.isna(sub.iloc[0]["sim"]):
                    sim_value = float(r["value"]) + 1000.0
                else:
                    sim_value = float(sub.iloc[0]["sim"])
                rows.append(dict(exp_id=exp, treatment=int(r["treatment"]),
                                 user_var=user_var, dssat=dssat_var, kind=r["kind"],
                                 date=pd.NaT, obs=float(r["value"]), sim=sim_value))
                seen_scalar.add(key)

    # time-series from PlantGro matched to FileT/CSV obs
    if obs_table is not None and not obs_table.empty:
        ts_obs = obs_table[obs_table["kind"] == "timeseries"]
        for exp, res in results.items():
            pg = getattr(res, "plantgro", None)
            o = ts_obs[ts_obs["exp_id"] == exp]
            if o.empty:
                continue
            oavg = o.groupby(["treatment", "date", "variable"], as_index=False)["value"].mean()
            for _, r in oavg.iterrows():
                col = r["variable"]
                if col not in ts_inv:
                    continue
                if pg is None or pg.empty or col not in pg.columns:
                    simv = np.nan
                else:
                    simv = _timeseries_sim_value(pg, int(r["treatment"]), r["date"], col)
                if pd.isna(simv):
                    simv = float(r["value"]) + float(
                        (cfg.get("objective", {}) or {}).get("missing_simulation_penalty", 1000.0)
                    )
                rows.append(dict(exp_id=exp, treatment=int(r["treatment"]),
                                 user_var=ts_inv.get(col, col), dssat=col, kind="timeseries",
                                 date=r["date"], obs=float(r["value"]), sim=float(simv)))

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = _drop_configured_zero_observations(df, cfg)
    if df.empty:
        return df

    # Attempt to join original sigmas/weights from obs_table if they exist
    if obs_table is not None and not obs_table.empty and "sigma" in obs_table.columns:
        lookup = {}
        for _, r in obs_table.iterrows():
            obs_var = str(r["variable"])
            if str(r.get("kind", "")) in ("scalar", "phenology"):
                _, obs_var = _scalar_observation_mapping(obs_var, sc_map, sc_inv)
            if pd.isna(r["date"]):
                key = (r["exp_id"], int(r["treatment"]), obs_var)
            else:
                key = (r["exp_id"], int(r["treatment"]), obs_var, pd.Timestamp(r["date"]).date())
            lookup[key] = (r["sigma"], r["weight"])
            
        def get_obs_params(row):
            var = row["dssat"]
            if pd.isna(row["date"]):
                key = (row["exp_id"], int(row["treatment"]), var)
            else:
                key = (row["exp_id"], int(row["treatment"]), var, pd.Timestamp(row["date"]).date())
                
            val = lookup.get(key)
            if val is not None:
                sig, wt = val
                if pd.isna(sig):
                    sig = _sigma(row["user_var"], row["obs"], cfg)
                if pd.isna(wt):
                    wts = cfg.get("objective", {}).get("weights", {}) or {}
                    wt = float(wts.get(row["user_var"], 1.0))
                return sig, wt
            
            sig = _sigma(row["user_var"], row["obs"], cfg)
            wts = cfg.get("objective", {}).get("weights", {}) or {}
            wt = float(wts.get(row["user_var"], 1.0))
            return sig, wt
            
        params = [get_obs_params(r) for _, r in df.iterrows()]
        df["sigma"] = [p[0] for p in params]
        df["weight"] = [p[1] for p in params]
    else:
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

    resid = resid.copy()
    resid["_loss"] = _weighted_loss(resid, cfg)
    loglik = float(-0.5 * np.sum(resid["_loss"] * resid["weight"]))

    if weighting == "sigma":
        # Total sigma-weighted chi-square (NO per-variable averaging). The
        # statistically 'correct' misfit; matches the Gaussian log-likelihood, and
        # is what the Bayesian engines expect.
        sc = float(np.sum(resid["_loss"] * resid["weight"]))
    elif weighting == "user":
        # Raw normalised RMSE per variable x explicit weights. Full manual control,
        # no automatic balancing of counts or scales.
        sc = 0.0
        for uv, g in resid.groupby("user_var"):
            sc += float(wts.get(uv, 1.0)) * _group_loss(g, cfg)
    elif weighting == "count_scale":
        # Each variable contributes the AVERAGE normalised squared error of its own
        # points (count-balanced: a 100-point series can't drown a 1-point yield),
        # then we average ACROSS variables -> a clean 'mean per-variable misfit'.
        per = [float(wts.get(uv, 1.0)) * _group_loss(g, cfg)
               for uv, g in resid.groupby("user_var")]
        sc = float(np.mean(per)) if per else float("inf")
    else:  # "unified" (default) and "agmip_wls"
        # SUM of per-variable mean normalised squared errors x weights. Count- and
        # scale-balanced like count_scale, but summed (so fitting more variables
        # costs more). 'agmip_wls' uses this same surface after the driver has reset
        # the per-variable weights to 1 / residual-variance (see pipeline).
        sc = 0.0
        for uv, g in resid.groupby("user_var"):
            sc += float(wts.get(uv, 1.0)) * _group_loss(g, cfg)

    penalty_cfg = (cfg.get("objective", {}) or {}).get("max_bias_penalty", {}) or {}
    if penalty_cfg:
        variable = str(penalty_cfg.get("variable", "anthesis"))
        lam = float(penalty_cfg.get("lambda", penalty_cfg.get("weight", 0.0)) or 0.0)
        if lam > 0:
            g = resid[resid["user_var"].astype(str) == variable]
            if not g.empty:
                max_abs = float(np.nanmax(np.abs(pd.to_numeric(g["resid"], errors="coerce"))))
                tolerance = float(penalty_cfg.get("tolerance", penalty_cfg.get("target", 0.0)) or 0.0)
                sigma = float(penalty_cfg.get("sigma", 1.0) or 1.0)
                excess = max(0.0, max_abs - tolerance)
                power = float(penalty_cfg.get("power", 2.0))
                sc += lam * float((excess / max(sigma, 1e-12)) ** power)

    per_var = {uv: metrics(g["obs"], g["sim"]) for uv, g in resid.groupby("user_var")}
    try:
        pev = (resid.groupby(["exp_id", "user_var"])
               .apply(lambda g: pd.Series(metrics(g["obs"], g["sim"])), include_groups=False)
               .reset_index())
    except TypeError:
        pev = (resid.groupby(["exp_id", "user_var"])
               .apply(lambda g: pd.Series(metrics(g["obs"], g["sim"])))
               .reset_index())
    return ObjectiveResult(score=float(sc), loglik=loglik, residuals=resid,
                           per_var=per_var, per_exp_var=pev)
