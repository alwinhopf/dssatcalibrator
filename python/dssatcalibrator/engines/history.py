"""Bayesian history matching for sparse cultivar/species calibration.

History matching is deliberately more cautious than "find the best theta".  It
asks which parameter sets are still plausible given the observations and their
uncertainty, then refocuses later waves inside that not-ruled-out-yet (NROY)
region.  This is a good fit when data are too sparse to identify a full
posterior, which is common for new DSSAT cultivars and analog species.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..samplers import sample


@dataclass
class HistoryResult:
    design: pd.DataFrame
    behavioural: pd.DataFrame
    best_theta: dict
    best_sample_id: int
    threshold: float
    ess: float
    obj_results: dict
    best: object
    waves: list = field(default_factory=list)


def _implausibility(result, mode: str = "max_z") -> float:
    resid = getattr(result, "residuals", None)
    if resid is None or resid.empty:
        return float("inf")
    if mode == "score":
        return float(np.sqrt(result.score)) if np.isfinite(result.score) else float("inf")
    z = np.abs(resid["resid"].to_numpy(dtype=float) / resid["sigma"].to_numpy(dtype=float))
    return float(np.nanmax(z)) if z.size else float("inf")


def _sample_box(space, n: int, low: np.ndarray, high: np.ndarray, seed: int) -> pd.DataFrame:
    unit = sample(space, n=n, engine="lhs", seed=seed, include_start=False)
    u = (unit.to_numpy(dtype=float) - space.low) / np.where(space.high > space.low, space.high - space.low, 1.0)
    vals = low + np.clip(u, 0.0, 1.0) * (high - low)
    return pd.DataFrame(vals, columns=space.names)


def run_history_matching(cfg: dict, score_results, space, *, progress: bool = True) -> HistoryResult:
    """Run iterative history matching and return the final NROY design."""
    hcfg = cfg.get("method", {}).get("bayesian", {}) or {}
    waves = int(hcfg.get("waves", 3))
    n = int(hcfg.get("n", hcfg.get("n_per_wave", 128)))
    cutoff = float(hcfg.get("implausibility_cutoff", 3.0))
    mode = str(hcfg.get("implausibility", "max_z")).lower()
    refocus_quantile = float(hcfg.get("refocus_quantile", 0.9))
    pad_frac = float(hcfg.get("pad_frac", 0.10))
    seed = int(cfg["calibrator"].get("seed", 42))

    low = np.asarray(space.low, dtype=float).copy()
    high = np.asarray(space.high, dtype=float).copy()
    all_rows, obj_results, wave_info = [], {}, []
    sid = 0

    for wave in range(waves):
        design = _sample_box(space, n, low, high, seed + wave)
        if wave == 0:
            start = pd.DataFrame([space.start], columns=space.names)
            design = pd.concat([start, design], ignore_index=True)
        thetas = [space.to_theta(design.iloc[i].to_numpy()) for i in range(len(design))]
        results = list(score_results(thetas))
        if not any(np.isfinite(r.score) for r in results):
            raise ValueError(
                "History matching found no valid candidates: every evaluated score "
                "is non-finite. Inspect the spawn manifest and per-run DSSAT errors."
            )
        impl = np.array([_implausibility(r, mode) for r in results], dtype=float)
        nroy_mask = impl <= cutoff
        if not nroy_mask.any() and np.isfinite(impl).any():
            finite = impl[np.isfinite(impl)]
            adaptive = float(np.quantile(finite, min(max(refocus_quantile, 0.0), 1.0)))
            nroy_mask = impl <= adaptive
        for i, (theta, res) in enumerate(zip(thetas, results)):
            obj_results[sid] = res
            all_rows.append({"sample_id": sid, "wave": wave, **theta,
                             "score": res.score, "loglik": res.loglik,
                             "n_obs": len(res.residuals),
                             "implausibility": float(impl[i]),
                             "nroy": bool(nroy_mask[i])})
            sid += 1
        nroy = design.loc[nroy_mask]
        wave_info.append({"wave": wave, "n": len(design), "nroy": int(nroy_mask.sum()),
                          "low": low.copy(), "high": high.copy()})
        if progress:
            print(f"  history wave {wave+1}/{waves}: NROY {int(nroy_mask.sum())}/{len(design)}", flush=True)
        if nroy.empty:
            break
        span = np.maximum(nroy.max().to_numpy() - nroy.min().to_numpy(), 1e-12)
        low = np.maximum(space.low, nroy.min().to_numpy() - pad_frac * span)
        high = np.minimum(space.high, nroy.max().to_numpy() + pad_frac * span)

    out = pd.DataFrame(all_rows)
    nroy_final = out[out["nroy"]].copy()
    if nroy_final.empty:
        valid = out[np.isfinite(out["score"])]
        nroy_final = valid.nsmallest(max(1, min(10, len(valid))), "score").copy()
    out["weight"] = 0.0
    if not nroy_final.empty:
        out.loc[nroy_final.index, "weight"] = 1.0 / len(nroy_final)
    valid = out[np.isfinite(out["score"])]
    if valid.empty:  # Defensive guard if a future wave policy changes.
        raise ValueError(
            "History matching found no valid candidates: every evaluated score "
            "is non-finite. Inspect the spawn manifest and per-run DSSAT errors."
        )
    best_sample_id = int(valid["score"].idxmin())
    best_theta = {name: float(out.loc[best_sample_id, name]) for name in space.names}
    best = obj_results[best_sample_id]
    return HistoryResult(design=out, behavioural=nroy_final, best_theta=best_theta,
                         best_sample_id=best_sample_id, threshold=cutoff,
                         ess=float(len(nroy_final)), obj_results=obj_results,
                         best=best, waves=wave_info)
