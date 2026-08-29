"""Surrogate / emulator acceleration (``method.surrogate``).

DSSAT is slow; a *surrogate* is a fast statistical stand-in. The idea:

1. Run a modest space-filling design on the **real** model (e.g. 64 runs).
2. Fit a cheap emulator ``score = f(parameters)`` (a Gaussian Process or Random
   Forest) to those runs.
3. Search the **emulator** over thousands of candidates for free, and keep the
   ``top_k`` most promising.
4. **Validate** only those ``top_k`` on the real DSSAT model.

You spend real model runs where they matter and let the emulator do the rest — this
is what makes the expensive MCMC/Sobol stages tractable for slow crops.

Returns the union of all *real* evaluations (training design + validated
candidates), which the orchestrator then post-processes with GLUE exactly like a
normal sampling run, so the figures and posterior summary are unchanged.

Requires scikit-learn (``pip install scikit-learn`` or ``dssatcalibrator[full]``).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..samplers import sample


@dataclass
class SurrogateResult:
    design: pd.DataFrame          # all REAL evaluations (train + validated)
    obj_results: dict             # sample_id -> ObjectiveResult
    info: dict = field(default_factory=dict)


def _fit_emulator(model: str, Xn: np.ndarray, y: np.ndarray, seed: int):
    if model == "rf":
        from sklearn.ensemble import RandomForestRegressor
        est = RandomForestRegressor(n_estimators=300, random_state=seed, n_jobs=-1)
    else:  # "gp"
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
        kernel = ConstantKernel(1.0) * RBF(length_scale=np.repeat(0.2, Xn.shape[1])) + WhiteKernel(1e-3)
        est = GaussianProcessRegressor(kernel=kernel, normalize_y=True,
                                       n_restarts_optimizer=2, random_state=seed)
    est.fit(Xn, y)
    return est


def run_surrogate(cfg: dict, space, score_results, *, progress: bool = True) -> SurrogateResult:
    """Emulator-accelerated calibration. ``score_results(list_of_theta) -> list of
    ObjectiveResult`` is the framework's parallel evaluator."""
    try:
        import sklearn  # noqa: F401
    except ImportError as exc:  # pragma: no cover - optional install
        raise ImportError(
            "The surrogate engine needs scikit-learn. Install it with "
            "`pip install scikit-learn` (or `pip install dssatcalibrator[full]`)."
        ) from exc

    scfg = cfg.get("method", {}).get("surrogate", {})
    model = str(scfg.get("engine", scfg.get("model", "gp"))).lower()
    n_train = int(scfg.get("n_train", 64))
    n_candidates = int(scfg.get("n_candidates", 5000))
    top_k = int(scfg.get("top_k", 10))
    seed = int(cfg["calibrator"].get("seed", 42))

    lo, hi = space.low, space.high
    span = np.where((hi - lo) == 0, 1.0, hi - lo)

    # 1-2. Train design on the real model.
    if progress:
        print(f"Surrogate: running {n_train} training points on the real model...", flush=True)
    train = sample(space, n=n_train, engine="lhs", seed=seed, include_start=True)
    train_thetas = [space.to_theta(train.iloc[i].to_numpy()) for i in range(len(train))]
    train_res = score_results(train_thetas)

    Xn, y, kept = [], [], []
    for th, r in zip(train_thetas, train_res):
        if np.isfinite(r.score):
            Xn.append([(th[n] - lo[j]) / span[j] for j, n in enumerate(space.names)])
            y.append(r.score)
            kept.append((th, r))
    if len(kept) < 4:
        raise RuntimeError("surrogate: too few successful training runs to fit an emulator")
    Xn = np.array(Xn); y = np.array(y)

    # 3. Fit emulator and score many candidates for free.
    if progress:
        print(f"Surrogate: fitting {model.upper()} emulator on {len(y)} runs, "
              f"searching {n_candidates} candidates...", flush=True)
    est = _fit_emulator(model, Xn, y, seed)
    cand = sample(space, n=n_candidates, engine="lhs", seed=seed + 1, include_start=False)
    cand_native = cand.to_numpy()
    cand_norm = (cand_native - lo) / span
    pred = est.predict(cand_norm)
    order = np.argsort(pred)[:top_k]

    # 4. Validate the most promising candidates on the real model.
    if progress:
        print(f"Surrogate: validating top {top_k} candidates on the real model...", flush=True)
    val_thetas = [space.to_theta(cand_native[i]) for i in order]
    val_res = score_results(val_thetas)

    # Combine all REAL evaluations into one design table.
    rows, obj_results = [], {}
    sid = 0
    for th, r in list(kept) + list(zip(val_thetas, val_res)):
        obj_results[sid] = r
        rows.append({"sample_id": sid, **th, "score": r.score, "loglik": r.loglik,
                     "n_obs": len(r.residuals)})
        sid += 1
    design = pd.DataFrame(rows)

    info = {"model": model, "n_train": int(len(y)), "n_validated": int(top_k),
            "best_predicted_score": float(pred[order[0]])}
    return SurrogateResult(design=design, obj_results=obj_results, info=info)
