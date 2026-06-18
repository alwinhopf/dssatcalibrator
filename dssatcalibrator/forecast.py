"""In-season LAI nowcast/forecast from a calibrated posterior (optional, OFF by default).

Once the model is calibrated (a posterior cloud of parameter sets, or a single best
fit), DSSAT already simulates the whole season — so a forward LAI estimate is just
the calibrated run read *past* the last observation. This module turns that into a
proper product:

1. **Ensemble** — run the behavioural parameter sets forward and take daily LAI
   percentiles (P10/P50/P90), so the forecast carries parameter uncertainty.
2. **Anchor continuity** — the simulated LAI on the last-observation date rarely
   equals the observed value; ``anchor_correction`` shifts the forward curve to
   start from the observation and decays the correction back toward the pure model
   over ``decay_days``. This gives a seam-free nowcast without needing to inject
   state into DSSAT (which CSM does not support cleanly).
3. **Skill** — ``lead_time_table`` summarises the ensemble spread vs lead time so a
   user knows how far past the last cloud-free image the estimate stays trustworthy.

Switch on with::

    forecast:
      active: true
      variables: ["LAID"]
      n_ensemble: 30           # behavioural sets to propagate (0 => best fit only)
      anchor_continuity: true
      decay_days: 21

The pure helpers (percentiles, anchor, skill) are unit-tested offline; the ensemble
propagation needs DSSAT runs (reuses ``orchestrator.spawn_results_for``).
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def ensemble_percentiles(curves: list[pd.DataFrame], variable: str = "LAID",
                         quantiles=(0.1, 0.5, 0.9)) -> pd.DataFrame:
    """Daily percentiles of ``variable`` across an ensemble of PlantGro curves.

    Each curve is a PlantGro DataFrame with ``date`` and ``variable`` columns.
    Returns a tidy frame: ``date | p10 | p50 | p90 | mean | n``.
    """
    frames = []
    for c in curves:
        if c is None or c.empty or variable not in c.columns or "date" not in c.columns:
            continue
        frames.append(c[["date", variable]].dropna())
    if not frames:
        return pd.DataFrame(columns=["date", "p10", "p50", "p90", "mean", "n"])
    allc = pd.concat(frames, ignore_index=True)
    g = allc.groupby("date")[variable]
    out = pd.DataFrame({
        "date": sorted(allc["date"].unique()),
    }).set_index("date")
    qs = g.quantile(list(quantiles)).unstack()
    out["p10"] = qs[quantiles[0]]
    out["p50"] = qs[quantiles[1]]
    out["p90"] = qs[quantiles[2]]
    out["mean"] = g.mean()
    out["n"] = g.count()
    return out.reset_index()


def anchor_correction(forecast: pd.DataFrame, last_obs_value: float, last_obs_date,
                      *, decay_days: int = 21, mode: str = "additive",
                      cols=("p10", "p50", "p90", "mean")) -> pd.DataFrame:
    """Shift a forecast to start from the last observation, decaying over ``decay_days``.

    At ``last_obs_date`` the correction fully closes the model–obs gap; it then
    relaxes linearly to zero by ``last_obs_date + decay_days``. ``additive`` shifts
    by the residual, ``multiplicative`` scales by the ratio. Pure function.
    """
    out = forecast.copy()
    if out.empty or pd.isna(last_obs_value):
        return out
    last_obs_date = pd.Timestamp(last_obs_date)
    anchor_row = out.loc[out["date"] == last_obs_date]
    if anchor_row.empty:
        # nearest prior date as anchor
        prior = out.loc[out["date"] <= last_obs_date]
        if prior.empty:
            return out
        anchor_row = prior.iloc[[-1]]
    sim_at_anchor = float(anchor_row["p50"].iloc[0])
    if mode == "multiplicative":
        if sim_at_anchor == 0:
            return out
        full = last_obs_value / sim_at_anchor
    else:
        full = last_obs_value - sim_at_anchor

    days_past = (out["date"] - last_obs_date).dt.days.clip(lower=0)
    weight = np.where(days_past <= decay_days, 1.0 - days_past / max(decay_days, 1), 0.0)
    weight = np.where(out["date"] < last_obs_date, 0.0, weight)  # only correct forward
    for c in cols:
        if c not in out.columns:
            continue
        if mode == "multiplicative":
            factor = 1.0 + (full - 1.0) * weight
            out[c + "_adj"] = out[c] * factor
        else:
            out[c + "_adj"] = out[c] + full * weight
    out["anchor_weight"] = weight
    return out


def lead_time_table(forecast: pd.DataFrame, last_obs_date) -> pd.DataFrame:
    """Forecast spread (P90−P10) and relative spread vs days past the last observation."""
    out = forecast.copy()
    if out.empty:
        return pd.DataFrame(columns=["lead_days", "p50", "spread", "rel_spread"])
    last_obs_date = pd.Timestamp(last_obs_date)
    out = out.loc[out["date"] >= last_obs_date].copy()
    out["lead_days"] = (out["date"] - last_obs_date).dt.days
    out["spread"] = out["p90"] - out["p10"]
    out["rel_spread"] = out["spread"] / out["p50"].replace(0, np.nan)
    return out[["lead_days", "date", "p50", "spread", "rel_spread"]].reset_index(drop=True)


def _behavioural_thetas(result, n: int) -> list[dict]:
    """Pick up to ``n`` behavioural parameter sets from a calibration result.

    Prefers the GLUE/SMC behavioural design (best-scoring rows); falls back to the
    single best theta. Returns a list of ``{param: value}`` dicts.
    """
    space = result.space
    design = getattr(result, "design", None)
    if n <= 0 or design is None or design.empty:
        return [dict(result.best_theta)]
    cols = [c for c in space.names if c in design.columns]
    if not cols:
        return [dict(result.best_theta)]
    ranked = design.sort_values("score") if "score" in design.columns else design
    thetas = [dict(zip(cols, ranked.iloc[i][cols])) for i in range(min(n, len(ranked)))]
    # always include the best fit
    if result.best_theta not in thetas:
        thetas.insert(0, dict(result.best_theta))
    return thetas[:max(n, 1)]


def forecast_lai(cfg: dict, result, *, last_obs: dict | None = None,
                 variable: str = "LAID") -> dict[str, pd.DataFrame]:
    """Propagate the calibrated ensemble forward and build per-experiment LAI forecasts.

    Parameters
    ----------
    last_obs
        optional ``{exp_id: (date, value)}`` of the last observation used for the
        anchor-continuity correction.

    Returns ``{exp_id: forecast_df}``. Requires DSSAT (runs the ensemble); the
    numeric post-processing is the unit-tested part above.
    """
    from .orchestrator import spawn_results_for

    fcfg = cfg.get("forecast", {})
    n_ens = int(fcfg.get("n_ensemble", 0))
    anchor = bool(fcfg.get("anchor_continuity", True))
    decay = int(fcfg.get("decay_days", 21))
    last_obs = last_obs or {}

    thetas = _behavioural_thetas(result, n_ens)
    logger.info("Forecast: propagating %d parameter set(s) for %s", len(thetas), variable)

    out: dict[str, pd.DataFrame] = {}
    for exp in result.experiments:
        curves = []
        for th in thetas:
            spawns = spawn_results_for(cfg, th, [exp])
            pg = spawns[exp].plantgro
            if not pg.empty:
                curves.append(pg)
        fc = ensemble_percentiles(curves, variable=variable)
        if anchor and exp in last_obs:
            d, v = last_obs[exp]
            fc = anchor_correction(fc, v, d, decay_days=decay)
        out[exp] = fc
    return out
