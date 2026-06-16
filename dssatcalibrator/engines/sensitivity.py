"""Sensitivity screening — *which parameters actually matter?* (``method.sensitivity``).

Before spending thousands of model runs calibrating 13 parameters, it pays to ask
which ones the output is even sensitive to. Screening ranks parameters by their
influence on the objective so you can **calibrate the few that matter and freeze
the rest** (use :func:`influential_params` + ``auto_activate``). This is the
recommended *first* stage for a new crop or site.

Methods
-------
``morris`` (default, cheap)
    *Elementary effects.* Walk through parameter space one parameter at a time and
    measure how much each step changes the output. Reported per parameter:

    * ``mu_star`` — average absolute effect = overall influence (rank by this);
    * ``mu``      — average signed effect (direction);
    * ``sigma``   — variability of the effect = interactions / non-linearity.

    Cost ≈ ``trajectories * (k + 1)`` model runs (``k`` = #parameters). Pure NumPy,
    no extra dependencies.

``sobol`` (variance-based, thorough, needs SALib)
    Decomposes output *variance* into each parameter's share. Reports first-order
    ``S1`` (effect alone) and total ``ST`` (effect including interactions). More
    rigorous, but costs many more runs. Requires ``pip install SALib`` (or the
    ``dssatcalibrator[full]`` extra).

All methods evaluate candidates through the framework's parallel evaluator, so the
whole screening design runs across every core in one batch.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import pandas as pd


@dataclass
class SensitivityResult:
    method: str
    ranking: pd.DataFrame        # one row per parameter, sorted most→least influential
    n_eval: int

    def influential(self, keep: int | None = None, rel_threshold: float = 0.1) -> list[str]:
        return influential_params(self.ranking, keep=keep, rel_threshold=rel_threshold)


def _to_native(unit_pts: np.ndarray, low: np.ndarray, high: np.ndarray) -> np.ndarray:
    """Map points from the unit hypercube [0,1]^k to native parameter units."""
    return low + unit_pts * (high - low)


def run_morris(space, score_results: Callable[[list[dict]], list], *,
               trajectories: int = 10, levels: int = 4, seed: int = 42,
               progress: bool = False) -> SensitivityResult:
    """Morris elementary-effects screening (pure NumPy)."""
    k = space.ndim
    rng = np.random.default_rng(seed)
    delta = levels / (2.0 * (levels - 1))            # standard Morris step
    grid = np.linspace(0.0, 1.0, levels)
    bases = grid[grid <= 1.0 - delta + 1e-9]         # bases that leave room for a +delta step

    # Build `trajectories` one-at-a-time walks of (k+1) points each.
    traj_pts, traj_order = [], []
    for _ in range(trajectories):
        x = rng.choice(bases, size=k)
        order = rng.permutation(k)
        pts = [x.copy()]
        cur = x.copy()
        for j in order:
            cur = cur.copy()
            cur[j] += delta
            pts.append(cur)
        traj_pts.append(np.array(pts))               # (k+1, k)
        traj_order.append(order)

    # Evaluate every point in one parallel batch.
    all_unit = np.vstack(traj_pts)                   # (trajectories*(k+1), k)
    thetas = [space.to_theta(row) for row in _to_native(all_unit, space.low, space.high)]
    results = score_results(thetas)
    scores = np.array([r.score if np.isfinite(r.score) else np.nan for r in results])

    # Elementary effects per parameter.
    ee = {name: [] for name in space.names}
    m = k + 1
    for t in range(trajectories):
        y = scores[t * m:(t + 1) * m]
        for step, j in enumerate(traj_order[t]):
            d = (y[step + 1] - y[step]) / delta
            if np.isfinite(d):
                ee[space.names[j]].append(d)

    rows = []
    for name in space.names:
        arr = np.array(ee[name], dtype=float)
        rows.append({"parameter": name,
                     "mu_star": float(np.mean(np.abs(arr))) if arr.size else np.nan,
                     "mu": float(np.mean(arr)) if arr.size else np.nan,
                     "sigma": float(np.std(arr)) if arr.size else np.nan})
    ranking = pd.DataFrame(rows).sort_values("mu_star", ascending=False).reset_index(drop=True)
    return SensitivityResult(method="morris", ranking=ranking, n_eval=len(thetas))


def run_sobol(space, score_results: Callable[[list[dict]], list], *,
              n_base: int = 256, seed: int = 42, progress: bool = False) -> SensitivityResult:
    """Sobol variance-based sensitivity (requires SALib)."""
    try:
        from SALib.analyze import sobol as sobol_analyze
        from SALib.sample import sobol as sobol_sample
    except ImportError as exc:  # pragma: no cover - depends on optional install
        raise ImportError(
            "Sobol sensitivity needs SALib. Install it with `pip install SALib` "
            "(or `pip install dssatcalibrator[full]`), or use method 'morris'."
        ) from exc

    problem = {"num_vars": space.ndim, "names": list(space.names),
               "bounds": list(zip(space.low.tolist(), space.high.tolist()))}
    X = sobol_sample.sample(problem, n_base)         # (n_base*(2k+2), k)
    thetas = [space.to_theta(row) for row in X]
    results = score_results(thetas)
    Y = np.array([r.score if np.isfinite(r.score) else np.nan for r in results])
    # SALib cannot handle NaN; replace failed runs with the worst observed score.
    if np.isnan(Y).any():
        Y = np.where(np.isnan(Y), np.nanmax(Y) if np.isfinite(np.nanmax(Y)) else 0.0, Y)
    Si = sobol_analyze.analyze(problem, Y, print_to_console=False)
    ranking = (pd.DataFrame({"parameter": problem["names"], "S1": Si["S1"], "ST": Si["ST"]})
               .sort_values("ST", ascending=False).reset_index(drop=True))
    return SensitivityResult(method="sobol", ranking=ranking, n_eval=len(thetas))


def run_sensitivity(space, score_results, *, method: str = "morris", **kwargs) -> SensitivityResult:
    """Dispatch to :func:`run_morris` or :func:`run_sobol`."""
    method = method.lower()
    if method == "morris":
        return run_morris(space, score_results,
                          trajectories=int(kwargs.get("trajectories", 10)),
                          levels=int(kwargs.get("levels", 4)),
                          seed=int(kwargs.get("seed", 42)),
                          progress=kwargs.get("progress", False))
    if method == "sobol":
        return run_sobol(space, score_results,
                         n_base=int(kwargs.get("n_base", 256)),
                         seed=int(kwargs.get("seed", 42)),
                         progress=kwargs.get("progress", False))
    raise ValueError(f"unknown sensitivity method '{method}' (use morris | sobol)")


def influential_params(ranking: pd.DataFrame, *, keep: int | None = None,
                       rel_threshold: float = 0.1) -> list[str]:
    """Pick the influential parameters from a ranking table.

    * If ``keep`` is given, return the top-``keep`` parameters.
    * Otherwise keep every parameter whose influence is at least ``rel_threshold``
      times the most influential one (default 10%).
    """
    metric = "mu_star" if "mu_star" in ranking.columns else "ST"
    r = ranking.sort_values(metric, ascending=False)
    if keep is not None:
        return r["parameter"].head(keep).tolist()
    top = float(r[metric].iloc[0]) if len(r) else 0.0
    if top <= 0:
        return r["parameter"].tolist()
    return r[r[metric] >= rel_threshold * top]["parameter"].tolist()


def anova_variance_share(design: pd.DataFrame, factor_cols: list[str],
                         response_col: str) -> pd.DataFrame:
    """Share of output variance attributable to each discrete factor (one-way ANOVA).

    Useful for the discrete-factor sweep (e.g. *which weather source best matches
    observations?*): given a table with one column per factor and a ``response_col``
    (e.g. the objective score), it returns each factor's between-group sum of
    squares as a fraction of the total — the classic "% of variance explained".
    """
    y = design[response_col].to_numpy(dtype=float)
    total_ss = float(np.sum((y - y.mean()) ** 2))
    rows = []
    for f in factor_cols:
        between = 0.0
        for _, g in design.groupby(f):
            gy = g[response_col].to_numpy(dtype=float)
            between += len(gy) * (gy.mean() - y.mean()) ** 2
        rows.append({"factor": f, "var_share": (between / total_ss) if total_ss > 0 else np.nan})
    return pd.DataFrame(rows).sort_values("var_share", ascending=False).reset_index(drop=True)
