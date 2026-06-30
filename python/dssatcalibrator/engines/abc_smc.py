"""ABC-SMC calibration for sparse or irregular observations.

Approximate Bayesian Computation uses the objective score as a distance between
simulation and observation.  It keeps particles below a tolerance and gradually
tightens that tolerance across waves.  It is useful when the likelihood is only a
rough convention, for example literature-derived species observations or mixed
phenology/yield summaries.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .. import priors
from ..samplers import sample


@dataclass
class AbcSmcResult:
    design: pd.DataFrame
    behavioural: pd.DataFrame
    best_theta: dict
    best_sample_id: int
    threshold: float
    ess: float
    obj_results: dict
    best: object
    thresholds: list = field(default_factory=list)
    initial_design: pd.DataFrame | None = None


def _theta_matrix(space, thetas: list[dict]) -> np.ndarray:
    return np.array([[t[n] for n in space.names] for t in thetas], dtype=float)


def _propose_from_particles(space, particles: list[dict], scores: np.ndarray,
                            n: int, rng: np.random.Generator) -> list[dict]:
    mat = _theta_matrix(space, particles)
    good = np.isfinite(scores)
    if good.sum() < 2:
        cov = np.diag(((space.high - space.low) * 0.05) ** 2)
    else:
        cov = np.cov(mat[good].T)
        if cov.ndim == 0:
            cov = np.array([[float(cov)]])
        cov = cov + np.diag(((space.high - space.low) * 0.01) ** 2)
    idx = rng.integers(0, len(particles), size=n)
    props = []
    for i in idx:
        base = np.array([particles[i][nm] for nm in space.names], dtype=float)
        draw = rng.multivariate_normal(base, cov)
        props.append(space.to_theta(np.clip(draw, space.low, space.high)))
    return props


def run_abc_smc(cfg: dict, score_results, space, *, progress: bool = True) -> AbcSmcResult:
    """Run a lightweight ABC-SMC population calibration."""
    bcfg = cfg.get("method", {}).get("bayesian", {}) or {}
    n_particles = int(bcfg.get("n_particles", 128))
    waves = int(bcfg.get("waves", 4))
    oversample = max(1, int(bcfg.get("oversample", 3)))
    quantile = float(bcfg.get("threshold_quantile", 0.5))
    min_accept = int(bcfg.get("min_accept", max(8, n_particles // 4)))
    seed = int(cfg["calibrator"].get("seed", 42))
    rng = np.random.default_rng(seed)

    if priors.has_informative_prior(space):
        init = priors.sample_prior_design(space, n_particles * oversample, rng)
    else:
        init = sample(space, n=n_particles * oversample, engine="lhs",
                      seed=seed, include_start=False)
    start = pd.DataFrame([space.start], columns=space.names)
    init = pd.concat([start, init], ignore_index=True)
    particles = [space.to_theta(init.iloc[i].to_numpy()) for i in range(len(init))]
    initial_design = pd.DataFrame([{"sample_id": i, **t} for i, t in enumerate(particles)])

    rows, obj_results, thresholds = [], {}, []
    sid = 0
    accepted_particles: list[dict] = []
    accepted_scores = np.array([], dtype=float)

    for wave in range(waves):
        if wave == 0:
            candidates = particles
        else:
            candidates = _propose_from_particles(space, accepted_particles, accepted_scores,
                                                 n_particles * oversample, rng)
        results = list(score_results(candidates))
        scores = np.array([r.score if np.isfinite(r.score) else np.inf for r in results], dtype=float)
        finite = scores[np.isfinite(scores)]
        if finite.size == 0:
            threshold = float("inf")
        elif wave == 0:
            threshold = float(np.quantile(finite, min(max(quantile, 0.0), 1.0)))
        else:
            previous = thresholds[-1]
            threshold = min(previous, float(np.quantile(finite, min(max(quantile, 0.0), 1.0))))
        keep = np.where(scores <= threshold)[0]
        if keep.size < min_accept and finite.size:
            keep = np.argsort(scores)[:min(len(scores), max(min_accept, n_particles))]
            threshold = float(scores[keep[-1]])
        keep = keep[np.argsort(scores[keep])[:n_particles]]
        accepted_particles = [candidates[i] for i in keep]
        accepted_scores = scores[keep]
        thresholds.append(threshold)

        for i, (theta, res) in enumerate(zip(candidates, results)):
            accepted = bool(i in set(keep.tolist()))
            obj_results[sid] = res
            rows.append({"sample_id": sid, "wave": wave, **theta,
                         "score": res.score, "loglik": res.loglik,
                         "n_obs": len(res.residuals), "threshold": threshold,
                         "accepted": accepted})
            sid += 1
        if progress:
            print(f"  ABC wave {wave+1}/{waves}: accepted {len(keep)}/{len(candidates)} "
                  f"(eps={threshold:.4g})", flush=True)
        if not accepted_particles:
            break

    design = pd.DataFrame(rows)
    design["weight"] = 0.0
    final = design[(design["wave"] == design["wave"].max()) & design["accepted"]].copy()
    if not final.empty:
        design.loc[final.index, "weight"] = 1.0 / len(final)
    valid = design[np.isfinite(design["score"])]
    best_sample_id = int(valid["score"].idxmin()) if not valid.empty else 0
    best_theta = {name: float(design.loc[best_sample_id, name]) for name in space.names}
    best = obj_results[best_sample_id]
    return AbcSmcResult(design=design, behavioural=final, best_theta=best_theta,
                        best_sample_id=best_sample_id,
                        threshold=float(thresholds[-1]) if thresholds else float("inf"),
                        ess=float(len(final)), obj_results=obj_results, best=best,
                        thresholds=thresholds, initial_design=initial_design)
