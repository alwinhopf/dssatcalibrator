"""Stepwise parameter selection — the AgMIP overfit guard (``method.select``).

The problem
-----------
Calibrating *every* parameter you can is a trap: with enough knobs you can fit the
calibration data beautifully and still predict new seasons badly (over-fitting,
a.k.a. *equifinality* — many parameter sets fit equally well). The AgMIP protocol
guards against this with a principled rule:

* **obligatory** parameters (``role: obligatory`` in the config) are *always*
  estimated — typically the degree-day / phenology coefficients.
* **candidate** parameters (``role: candidate``, the default for anything not marked
  obligatory) are added **one at a time, and kept only if doing so lowers an
  information criterion** that penalises extra parameters.

Information criteria (lower is better)::

    BIC  = k*ln(n) - 2*logL          # heavier penalty, prefers simpler models
    AICc = 2k - 2*logL + 2k(k+1)/(n-k-1)   # small-sample-corrected AIC

where ``k`` = number of estimated parameters, ``n`` = number of observations, and
``logL`` is the Gaussian log-likelihood at the best fit. A candidate is only added
if it reduces the criterion by more than ``min_delta``.

Cost note
---------
Each tentative parameter set is *fully optimised* (Nelder-Mead), so selection runs
many optimisations. Keep ``optimizer_restarts`` and ``maxiter`` modest, and treat
this as a deliberate, heavier analysis. The returned subset can then be fed to a
single richer run (SMC-PF / MCMC).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from ..spaces import ParameterSpace
from .optimizers import run_optimizer


@dataclass
class SelectionResult:
    selected: list[str]          # parameter names chosen for calibration
    obligatory: list[str]
    criterion: str               # "bic" | "aicc"
    best_theta: dict             # best fit over the selected set (others at start)
    history: list = field(default_factory=list)   # [{"step","added","value","k","n"}]


def _criterion(name: str, loglik: float, k: int, n: int) -> float:
    if not np.isfinite(loglik):
        return float("inf")
    if name == "aicc":
        aic = 2 * k - 2 * loglik
        denom = max(n - k - 1, 1)
        return aic + (2 * k * (k + 1)) / denom
    return k * np.log(max(n, 1)) - 2 * loglik   # bic (default)


def _subspace(space_full: ParameterSpace, names: list[str]) -> ParameterSpace:
    """A ParameterSpace restricted to a subset of parameter names (same bounds)."""
    keep = [s for s in space_full.specs if s["name"] in set(names)]
    return ParameterSpace(
        names=[s["name"] for s in keep],
        low=np.array([float(s["min"]) for s in keep]),
        high=np.array([float(s["max"]) for s in keep]),
        start=np.array([float(s.get("start", 0.5 * (float(s["min"]) + float(s["max"])))) for s in keep]),
        specs=keep,
    )


def stepwise_select(space_full: ParameterSpace,
                    score_results: Callable[[list[dict]], list],
                    *, criterion: str = "bic", optimizer: str = "nelder_mead",
                    optimizer_restarts: int = 2, maxiter: int | None = None,
                    min_delta: float = 0.0, seed: int = 42,
                    progress: bool = False) -> SelectionResult:
    """Run AgMIP stepwise selection.

    Parameters
    ----------
    space_full
        the full :class:`ParameterSpace` over *all* active parameters; each spec may
        carry ``role: obligatory | candidate``.
    score_results
        ``score_results(list_of_full_theta) -> list_of_ObjectiveResult`` — evaluates
        complete parameter vectors (held parameters at their ``start``) through the
        framework's parallel evaluator. Each result exposes ``.score``, ``.loglik``
        and ``.residuals`` (whose length is ``n``).
    """
    criterion = criterion.lower()
    start_full = {s["name"]: float(s.get("start", 0.5 * (float(s["min"]) + float(s["max"]))))
                  for s in space_full.specs}

    obligatory = [s["name"] for s in space_full.specs if str(s.get("role", "candidate")) == "obligatory"]
    candidates = [s["name"] for s in space_full.specs if s["name"] not in obligatory]
    # If nobody is marked obligatory, fall back to the first parameter so the search
    # always estimates at least one thing.
    if not obligatory:
        obligatory = [space_full.names[0]]
        candidates = [n for n in space_full.names if n != obligatory[0]]

    def _fit(names: list[str]):
        """Optimise over ``names`` (others held at start); return (best_full_theta, crit)."""
        held = {n: start_full[n] for n in space_full.names if n not in names}
        sub = _subspace(space_full, names)

        def score_batch(subset_thetas):
            full = [{**held, **st} for st in subset_thetas]
            return [r.score for r in score_results(full)]

        opt = run_optimizer(sub, score_batch, method=optimizer, seed=seed,
                            restarts=optimizer_restarts, maxiter=maxiter)
        best_full = {**held, **opt.best_theta}
        res = score_results([best_full])[0]
        n = len(res.residuals)
        crit = _criterion(criterion, res.loglik, k=len(names), n=n)
        return best_full, crit, n

    selected = list(obligatory)
    best_full, best_crit, n_obs = _fit(selected)
    history = [{"step": 0, "added": "+".join(obligatory), "value": round(best_crit, 3),
                "k": len(selected), "n": n_obs}]
    if progress:
        print(f"[select] obligatory {selected} -> {criterion.upper()}={best_crit:.3f}", flush=True)

    remaining = list(candidates)
    step = 0
    while remaining:
        step += 1
        trials = []
        for cand in remaining:
            theta_c, crit_c, n_c = _fit(selected + [cand])
            trials.append((cand, crit_c, theta_c))
            if progress:
                print(f"[select] try +{cand}: {criterion.upper()}={crit_c:.3f} "
                      f"(current {best_crit:.3f})", flush=True)
        cand, crit_c, theta_c = min(trials, key=lambda t: t[1])
        if crit_c < best_crit - min_delta:
            selected.append(cand)
            remaining.remove(cand)
            best_crit, best_full = crit_c, theta_c
            history.append({"step": step, "added": cand, "value": round(crit_c, 3),
                            "k": len(selected), "n": n_obs})
            if progress:
                print(f"[select] KEEP +{cand} -> {criterion.upper()}={crit_c:.3f}", flush=True)
        else:
            if progress:
                print(f"[select] stop: no candidate lowers {criterion.upper()}", flush=True)
            break

    return SelectionResult(selected=selected, obligatory=obligatory, criterion=criterion,
                           best_theta=best_full, history=history)
