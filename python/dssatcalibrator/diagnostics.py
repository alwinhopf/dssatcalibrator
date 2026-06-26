"""Identifiability and structural-adequacy diagnostics (optional, OFF by default).

Sparse-data calibration (a new cultivar/species from a few literature site-years)
is prone to two failure modes the headline fit metrics hide:

* **Non-identifiability / equifinality** — a parameter the data cannot constrain,
  whose posterior stays as wide as its prior, or that trades off with another. The
  fit looks fine but the coefficient is meaningless and won't transfer.
* **Structural inadequacy** — *no* parameter set can reproduce the observations
  (wrong analog module, missing process). Calibration then bends coefficients into
  unphysical values to compensate.

``identifiability`` and ``structural_adequacy`` surface both from an existing
calibration result. Both are pure (no DSSAT) and unit-tested.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .config import active_parameters

logger = logging.getLogger(__name__)


def _prior_std(spec: dict) -> float:
    """Approximate prior standard deviation for a parameter spec."""
    prior = spec.get("prior") or {}
    dist = str(prior.get("dist", "uniform")).lower()
    lo, hi = float(spec["min"]), float(spec["max"])
    if dist == "normal" and "sd" in prior:
        return float(prior["sd"])
    if dist == "lognormal" and "sd" in prior:
        return float(prior["sd"])
    return (hi - lo) / np.sqrt(12.0)            # uniform


def identifiability(result, *, behavioural_quantile: float = 0.1) -> pd.DataFrame:
    """Per-parameter identifiability from the behavioural (best-scoring) design.

    Returns a table with the posterior std, the prior std, their ratio (small =
    well identified), and the strongest pairwise correlation with another active
    parameter (|r| near 1 = the two are not separately identifiable).
    """
    cfg = result.cfg
    space = result.space
    design = getattr(result, "design", None)
    specs = {s["name"]: s for s in active_parameters(cfg)}
    names = [n for n in space.names if n in (design.columns if design is not None else [])]
    if design is None or design.empty or not names:
        return pd.DataFrame(columns=["parameter", "posterior_sd", "prior_sd",
                                     "sd_ratio", "max_abs_corr", "identifiable"])

    d = design
    if "score" in d.columns and len(d) > 5:
        k = max(5, int(np.ceil(behavioural_quantile * len(d))))
        d = d.nsmallest(k, "score")
    sub = d[names].astype(float)

    corr = sub.corr().abs()
    # ``DataFrame.values`` can be a read-only view (newer pandas/numpy), so write
    # the NaN diagonal on an owned copy and rebuild the frame.
    _carr = corr.to_numpy(copy=True)
    np.fill_diagonal(_carr, np.nan)
    corr = pd.DataFrame(_carr, index=corr.index, columns=corr.columns)

    rows = []
    for n in names:
        post_sd = float(sub[n].std(ddof=1)) if len(sub) > 1 else np.nan
        pri_sd = _prior_std(specs[n]) if n in specs else np.nan
        ratio = post_sd / pri_sd if (pri_sd and pri_sd > 0) else np.nan
        max_corr = float(np.nanmax(corr[n].values)) if n in corr else np.nan
        rows.append({
            "parameter": n,
            "posterior_sd": post_sd,
            "prior_sd": pri_sd,
            "sd_ratio": ratio,
            "max_abs_corr": max_corr,
            # heuristic: identifiable if posterior is meaningfully tighter than prior
            "identifiable": bool(np.isfinite(ratio) and ratio < 0.6),
        })
    out = pd.DataFrame(rows).sort_values("sd_ratio").reset_index(drop=True)

    weak = out.loc[~out["identifiable"], "parameter"].tolist()
    collinear = out.loc[out["max_abs_corr"] > 0.9, "parameter"].tolist()
    if weak:
        logger.warning("Weakly identified parameter(s) (posterior ~ prior): %s — "
                       "consider fixing them or adding contrasting data.", weak)
    if collinear:
        logger.warning("Strongly correlated (|r|>0.9) parameter(s): %s — not "
                       "separately identifiable; consider fixing one.", collinear)
    return out


def structural_adequacy(result, *, ef_floor: float = 0.0,
                        nrmse_ceiling: float = 50.0) -> pd.DataFrame:
    """Per-variable check that the *best* fit is structurally capable.

    Flags variables where the best-fit modelling efficiency EF is below ``ef_floor``
    (model no better than the observed mean) or nRMSE% exceeds ``nrmse_ceiling`` —
    signs the chosen module/analog cannot reproduce the data, which calibration
    cannot fix. Returns one row per variable with a ``flag`` and reason.
    """
    rows = []
    per_var = getattr(result.best, "per_var", {}) or {}
    for var, m in per_var.items():
        ef = m.get("EF", np.nan)
        nrmse = m.get("nRMSE_pct", np.nan)
        reasons = []
        if np.isfinite(ef) and ef < ef_floor:
            reasons.append(f"EF={ef:.2f}<{ef_floor}")
        if np.isfinite(nrmse) and nrmse > nrmse_ceiling:
            reasons.append(f"nRMSE={nrmse:.0f}%>{nrmse_ceiling}%")
        rows.append({"variable": var, "EF": ef, "nRMSE_pct": nrmse,
                     "n": m.get("n", np.nan),
                     "flag": bool(reasons), "reason": "; ".join(reasons)})
    out = pd.DataFrame(rows)
    flagged = out.loc[out["flag"], "variable"].tolist() if not out.empty else []
    if flagged:
        logger.warning("Possible structural inadequacy for %s — the best fit still "
                       "misses badly; check the analog module/observations before "
                       "trusting the calibration.", flagged)
    return out
