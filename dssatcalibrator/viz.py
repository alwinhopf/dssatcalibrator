"""Visualisation & reporting for a calibration run.

``make_report(result, outdir)`` writes a set of PNGs + CSVs:

* ``fig_param_posteriors`` — prior (sampled) vs posterior-weighted distribution
  per parameter, with start and best-fit marked.
* ``fig_score_funnel``     — score distribution + the spawn funnel
  (prior samples -> behavioural -> best) and posterior ESS.
* ``fig_obs_vs_sim``       — 1:1 simulated-vs-observed for the best fit.
* ``fig_timeseries``       — best-fit simulated growth curves vs observed points.
* ``fig_fit_bars``         — per-variable nRMSE / Willmott d.
* ``summary_fit.csv`` / ``posterior_summary.csv`` / ``design.csv`` / ``best_theta.json``.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_ACCENT = "#534AB7"
_OBS = "#D85A30"


def _grid(n):
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))
    return rows, cols


def plot_param_posteriors(result, path):
    space, design = result.space, result.design
    names = space.names
    w = design["weight"].to_numpy() if "weight" in design else np.ones(len(design))
    # For SMC the prior is the initial ensemble (design holds the moved/posterior
    # cloud); for GLUE the sampled design itself *is* the prior.
    prior = (result.extras or {}).get("initial_design")
    prior = prior if (prior is not None and not prior.empty) else design
    rows, cols = _grid(len(names))
    fig, axes = plt.subplots(rows, cols, figsize=(3.3 * cols, 2.6 * rows), squeeze=False)
    for k, name in enumerate(names):
        ax = axes[k // cols][k % cols]
        x = design[name].to_numpy(dtype=float)
        xp = prior[name].to_numpy(dtype=float)
        lo, hi = space.low[k], space.high[k]
        bins = np.linspace(lo, hi, 16)
        ax.hist(xp, bins=bins, density=True, color="#C9C7E8", alpha=0.8, label="prior")
        if w.sum() > 0:
            ax.hist(x, bins=bins, weights=w, density=True, histtype="step",
                    color=_ACCENT, linewidth=2, label="posterior")
        ax.axvline(space.start[k], color="#888780", ls=":", lw=1.5, label="start")
        ax.axvline(result.best_theta[name], color=_OBS, lw=2, label="best")
        ax.set_title(name, fontsize=10)
        ax.set_yticks([])
    for j in range(len(names), rows * cols):
        axes[j // cols][j % cols].axis("off")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.suptitle("Parameter prior vs posterior", y=1.02, fontsize=12)
    fig.tight_layout(rect=(0, 0.06, 1, 1.0))
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=9, frameon=False)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_score_funnel(result, path):
    design = result.design
    glue = result.glue
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    scores = design["score"].replace([np.inf, -np.inf], np.nan).dropna()
    ax1.hist(scores, bins=30, color="#C9C7E8", edgecolor="white")
    if glue is not None and np.isfinite(glue.threshold):
        ax1.axvline(glue.threshold, color=_OBS, lw=2,
                    label=f"behavioural threshold\n(q={result.cfg['method'].get('bayesian',{}).get('behavioural_quantile',0.1)})")
        ax1.legend(fontsize=9, frameon=False)
    ax1.set_xlabel("objective score (lower = better fit)")
    ax1.set_ylabel("spawn count")
    ax1.set_title("Score distribution across spawns")

    n_total = len(design)
    n_behav = len(glue.behavioural) if glue is not None else 0
    stages = ["prior\nsamples", "behavioural", "best"]
    counts = [n_total, n_behav, 1]
    ax2.bar(stages, counts, color=["#C9C7E8", _ACCENT, _OBS])
    for i, c in enumerate(counts):
        ax2.text(i, c, str(c), ha="center", va="bottom", fontsize=10)
    ess = getattr(glue, "ess", float("nan"))
    ax2.set_title(f"Spawn funnel  (posterior ESS = {ess:.1f})")
    ax2.set_ylabel("number of spawns")
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_ess_trajectory(result, path):
    """SMC only: effective sample size over the sequential assimilation steps,
    with the resample/move threshold and the steps where resampling fired."""
    trace = getattr(result.glue, "ess_trace", None)
    if not trace:
        return
    df = pd.DataFrame(trace)
    n = int(df["n"].iloc[0])
    ess_frac = float(result.cfg.get("method", {}).get("bayesian", {}).get("ess_frac", 0.5))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(df["step"], df["ess"], marker="o", color=_ACCENT, lw=1.8, label="ESS")
    ax.axhline(n * ess_frac, color=_OBS, ls="--", lw=1.5,
               label=f"resample threshold ({ess_frac:g}·N)")
    ax.axhline(n, color="#888780", ls=":", lw=1, label=f"N = {n}")
    fired = df[df["resampled"]]
    if not fired.empty:
        ax.scatter(fired["step"], fired["ess"], s=90, facecolors="none",
                   edgecolors=_OBS, linewidths=1.8, zorder=5, label="resample + move")
    ax.set_xlabel("assimilation step (time-series dates, then end-of-season scalars)")
    ax.set_ylabel("effective sample size")
    ax.set_ylim(0, n * 1.05)
    ax.set_title("SMC particle filter — ESS trajectory")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_obs_vs_sim(result, path):
    resid = result.best.residuals
    if resid.empty:
        return
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    vars_ = sorted(resid["user_var"].unique())
    cmap = plt.get_cmap("tab10")
    for i, uv in enumerate(vars_):
        g = resid[resid["user_var"] == uv]
        ax.scatter(g["obs"], g["sim"], s=28, alpha=0.8, color=cmap(i % 10), label=uv)
    lim = [0, max(resid["obs"].max(), resid["sim"].max()) * 1.05]
    ax.plot(lim, lim, ls="--", color="#888780", lw=1)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("observed"); ax.set_ylabel("simulated (best fit)")
    ax.set_title("Observed vs simulated — best fit")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_timeseries(result, best_spawns, path):
    cfg = result.cfg
    ts_map = cfg.get("engine", {}).get("timeseries_outputs", {})
    obs = result.obs.table
    obs_ts = obs[obs["kind"] == "timeseries"] if not obs.empty else obs
    if not ts_map or obs_ts.empty:
        return
    user_vars = list(ts_map.keys())
    rows, cols = _grid(len(user_vars))
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.2 * rows), squeeze=False)
    cmap = plt.get_cmap("tab10")
    for k, uv in enumerate(user_vars):
        ax = axes[k // cols][k % cols]
        col = ts_map[uv]
        for ei, (exp, res) in enumerate(best_spawns.items()):
            pg = getattr(res, "plantgro", None)
            if pg is None or pg.empty or col not in pg.columns:
                continue
            for trt, gpg in pg.groupby("treatment"):
                ax.plot(gpg["date"], gpg[col], color=cmap(ei % 10), lw=1, alpha=0.7)
            o = obs_ts[(obs_ts["exp_id"] == exp) & (obs_ts["variable"] == col)]
            if not o.empty:
                ax.scatter(o["date"], o["value"], color=cmap(ei % 10), s=22,
                           edgecolor="black", linewidth=0.4, label=exp, zorder=5)
        ax.set_title(f"{uv}  ({col})", fontsize=10)
        ax.tick_params(axis="x", labelrotation=30, labelsize=7)
    for j in range(len(user_vars), rows * cols):
        axes[j // cols][j % cols].axis("off")
    handles, labels = [], []
    for axrow in axes:
        for ax in axrow:
            h, l = ax.get_legend_handles_labels()
            for hh, ll in zip(h, l):
                if ll not in labels:
                    handles.append(hh); labels.append(ll)
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=min(6, len(labels)),
                   fontsize=8, frameon=False)
    fig.suptitle("Best-fit growth curves (lines) vs observations (points)", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_fit_bars(result, path):
    pv = result.best.per_var
    if not pv:
        return
    df = pd.DataFrame(pv).T.reset_index().rename(columns={"index": "variable"})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.bar(df["variable"], df["nRMSE_pct"], color=_ACCENT)
    ax1.set_ylabel("nRMSE (%)"); ax1.set_title("Relative error per variable")
    ax1.tick_params(axis="x", labelrotation=30)
    ax2.bar(df["variable"], df["d"], color="#1D9E75")
    ax2.set_ylabel("Willmott d"); ax2.set_ylim(0, 1); ax2.set_title("Agreement index per variable")
    ax2.tick_params(axis="x", labelrotation=30)
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def summary_fit_table(result) -> pd.DataFrame:
    pv = result.best.per_var
    if not pv:
        return pd.DataFrame()
    df = pd.DataFrame(pv).T.reset_index().rename(columns={"index": "variable"})
    cols = ["variable", "n", "RMSE", "nRMSE_pct", "MBE", "d", "EF", "R2"]
    return df[[c for c in cols if c in df.columns]].round(3)


def make_report(result, outdir, best_spawns=None, figdir=None) -> dict:
    """Write a calibration report. Data tables -> ``outdir``; all PNG figures ->
    ``figdir`` (defaults to ``outdir`` if not given). Returns the paths written."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir = Path(figdir) if figdir is not None else outdir
    figdir.mkdir(parents=True, exist_ok=True)
    paths = {}

    result.design.to_csv(outdir / "design.csv", index=False)
    paths["design"] = outdir / "design.csv"

    fit = summary_fit_table(result)
    fit.to_csv(outdir / "summary_fit.csv", index=False)
    paths["summary_fit"] = outdir / "summary_fit.csv"

    (outdir / "best_theta.json").write_text(json.dumps(result.best_theta, indent=2))
    if result.glue is not None:
        from .engines.glue import posterior_summary
        posterior_summary(result.glue, result.space.names).round(4).to_csv(
            outdir / "posterior_summary.csv", index=False)
        paths["posterior_summary"] = outdir / "posterior_summary.csv"

    if not result.best.per_exp_var.empty:
        result.best.per_exp_var.round(3).to_csv(outdir / "fit_by_experiment.csv", index=False)

    plot_param_posteriors(result, figdir / "fig_param_posteriors.png")
    plot_score_funnel(result, figdir / "fig_score_funnel.png")
    plot_ess_trajectory(result, figdir / "fig_ess_trajectory.png")
    plot_obs_vs_sim(result, figdir / "fig_obs_vs_sim.png")
    plot_fit_bars(result, figdir / "fig_fit_bars.png")
    if best_spawns:
        plot_timeseries(result, best_spawns, figdir / "fig_timeseries.png")

    if result.nsga2 is not None:
        _plot_pareto(result.nsga2, figdir / "fig_pareto.png")
        result.nsga2.front().round(4).to_csv(outdir / "pareto_front.csv", index=False)

    for key in ("param_posteriors", "score_funnel", "ess_trajectory", "obs_vs_sim",
                "fit_bars", "timeseries", "pareto"):
        p = figdir / f"fig_{key}.png"
        if p.exists():
            paths[key] = p
    return paths


def _plot_pareto(nsga2, path):
    F = nsga2.F
    vars_ = nsga2.objective_vars
    if F.shape[1] == 2:
        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        ax.scatter(F[:, 0], F[:, 1], color=_ACCENT, s=30)
        ax.set_xlabel(f"nRMSE {vars_[0]} (%)"); ax.set_ylabel(f"nRMSE {vars_[1]} (%)")
        ax.set_title("NSGA-II Pareto front")
    else:
        fig, ax = plt.subplots(figsize=(6.5, 4))
        for i in range(F.shape[0]):
            ax.plot(range(len(vars_)), F[i], color=_ACCENT, alpha=0.5, marker="o")
        ax.set_xticks(range(len(vars_))); ax.set_xticklabels(vars_, rotation=30)
        ax.set_ylabel("nRMSE (%)"); ax.set_title("NSGA-II Pareto front (parallel coords)")
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
