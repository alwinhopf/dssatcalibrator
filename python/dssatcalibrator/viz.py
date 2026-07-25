"""Visualisation & reporting for a calibration run.

``make_report(result, outdir)`` writes a set of PNGs + CSVs:

* ``fig_param_posteriors`` — prior (sampled) vs posterior-weighted distribution
  per parameter, with start and best-fit marked.
* ``fig_score_funnel``     — score distribution + the spawn funnel
  (prior samples -> behavioural -> best) and posterior ESS.
* ``fig_obs_vs_sim``       — 1:1 simulated-vs-observed for the best fit.
* ``fig_timeseries``       — best-fit simulated growth curves vs observed points.
* ``fig_fit_bars``         — per-variable nRMSE / Willmott d.
* ``summary_fit.csv`` / ``phenology_report.csv`` / ``objective_breakdown.csv`` /
  ``manifest.csv`` / ``posterior_summary.csv`` / ``design.csv`` /
  ``best_theta.json``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .config import crop_for, resolve_dssat_paths
from .spawn import theta_hash
from .weather import read_wth
from .writers import parse_fields

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


def plot_obs_vs_sim_by_category(result, path):
    resid = result.best.residuals
    if resid.empty:
        return
    user_vars = sorted(resid["user_var"].unique())
    n_vars = len(user_vars)
    if n_vars == 0:
        return
    rows, cols = _grid(n_vars)
    fig, axes = plt.subplots(rows, cols, figsize=(4.0 * cols, 4.0 * rows), squeeze=False)
    cmap = plt.get_cmap("tab10")
    for i, uv in enumerate(user_vars):
        ax = axes[i // cols][i % cols]
        g = resid[resid["user_var"] == uv]
        ax.scatter(g["obs"], g["sim"], s=30, alpha=0.8, color=cmap(i % 10), edgecolors="none")
        min_val = min(g["obs"].min(), g["sim"].min())
        max_val = max(g["obs"].max(), g["sim"].max())
        span = max_val - min_val
        if span == 0:
            span = 1.0
        lim = [min_val - 0.05 * span, max_val + 0.05 * span]
        ax.plot(lim, lim, ls="--", color="#888780", lw=1)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel("observed")
        ax.set_ylabel("simulated (best fit)")
        ax.set_title(uv, fontsize=11, fontweight="bold")
    for j in range(n_vars, rows * cols):
        axes[j // cols][j % cols].axis("off")
    fig.suptitle("Observed vs Simulated by Category", y=1.02, fontsize=13, fontweight="bold")
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


_EXPERIMENT_PANEL_SPECS = [
    {
        "title": "Leaf area index",
        "series": [("LAID", "LAI")],
    },
    {
        "title": "Growth stages",
        "stage": True,
    },
    {
        "title": "Canopy height and width",
        "series": [("CHTD", "height"), ("CWID", "width")],
    },
    {
        "title": "Total aboveground biomass",
        "series": [("CWAD", "aboveground")],
    },
    {
        "title": "Organ biomass",
        "series": [("SWAD", "stem"), ("LWAD", "leaf"), ("RWAD", "root"), ("GWAD", "grain")],
    },
    {
        "title": "Water and nitrogen stress",
        "series": [
            ("WSPD", "water photo"), ("WSDD", "water dev"), ("WSGD", "water growth"),
            ("NSTD", "N stress"), ("NSPD", "N photo"), ("NSGD", "N growth"),
        ],
    },
    {
        "title": "Soil water and nitrogen",
        "soil": True,
    },
    {
        "title": "Weather",
        "series": [("TMIN", "tmin"), ("TMAX", "tmax"), ("SRAD", "solar radiation")],
        "sources": ["Weather.OUT", "WTH"],
        "daily_average": True,
    },
    {
        "title": "Tissue nitrogen concentration",
        "series": [("LN%D", "leaf N%"), ("SN%D", "stem N%"), ("RN%D", "root N%"), ("GN%D", "grain N%")],
        "sources": ["PlantN.OUT", "PlantGro.OUT"],
    },
]


def _obs_aliases(result) -> dict[str, set[str]]:
    aliases: dict[str, set[str]] = {}
    engine = result.cfg.get("engine", {}) if hasattr(result, "cfg") else {}
    for block in ("timeseries_outputs", "scalar_outputs"):
        for user_var, dssat_var in (engine.get(block, {}) or {}).items():
            aliases.setdefault(str(dssat_var), set()).add(str(user_var))
            aliases.setdefault(str(user_var), set()).add(str(dssat_var))
    return aliases


def _obs_for_variable(obs: pd.DataFrame, exp_id: str, variable: str,
                      aliases: dict[str, set[str]], treatment: int | None = None) -> pd.DataFrame:
    if obs is None or obs.empty:
        return pd.DataFrame()
    wanted = {variable, *aliases.get(variable, set())}
    out = obs[(obs["exp_id"].astype(str) == str(exp_id)) & (obs["variable"].astype(str).isin(wanted))].copy()
    if treatment is not None and "treatment" in out.columns:
        out = out[pd.to_numeric(out["treatment"], errors="coerce") == int(treatment)]
    if "date" in out.columns:
        out["date"] = pd.to_datetime(out["date"], errors="coerce")
    return out


def _combined_output_long(spawn) -> pd.DataFrame:
    outputs = getattr(spawn, "outputs", None) or {}
    frames = []
    long = outputs.get("long") if isinstance(outputs, dict) else None
    if long is not None and not long.empty:
        frames.append(long.copy())
    pg = getattr(spawn, "plantgro", pd.DataFrame())
    if pg is not None and not pg.empty:
        value_cols = [
            c for c in pg.columns
            if c not in {"run", "treatment", "date"} and pd.api.types.is_numeric_dtype(pg[c])
        ]
        if value_cols:
            id_vars = [c for c in ("run", "treatment", "date", "DAP") if c in pg.columns]
            pgl = pg.melt(id_vars=id_vars, value_vars=value_cols,
                          var_name="variable", value_name="value")
            pgl["source_file"] = "PlantGro.OUT"
            frames.append(pgl)
    if not frames:
        return pd.DataFrame(columns=["source_file", "treatment", "date", "DAP", "variable", "value"])
    out = pd.concat(frames, ignore_index=True, sort=False)
    if "date" in out.columns:
        out["date"] = pd.to_datetime(out["date"], errors="coerce")
    if "treatment" in out.columns:
        out["treatment"] = pd.to_numeric(out["treatment"], errors="coerce").astype("Int64")
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    return out.dropna(subset=["value"])


def _filex_override_value(cfg: dict, exp_id: str, section: str, field: str):
    for rec in (cfg.get("filex_overrides", {}) or {}).get(str(exp_id), []) or []:
        if not isinstance(rec, dict):
            continue
        if str(rec.get("section", "")).upper() == section.upper() and str(rec.get("field", "")).upper() == field.upper():
            return rec.get("value")
    return None


def _weather_station_for(result, exp_id: str) -> str | None:
    cfg = getattr(result, "cfg", {}) or {}
    override = _filex_override_value(cfg, exp_id, "FIELDS", "WSTA")
    if override:
        return str(override)
    try:
        crop = crop_for(cfg, (cfg.get("crops") or [{}])[0].get("code", "HM"))
        filex = Path(cfg["source"]["hemp_dir"]) / f"{exp_id}.{crop['filex_ext']}"
        return parse_fields(filex).get("wsta")
    except Exception:
        return None


def _weather_fallback_long(result, exp_id: str, sim_long: pd.DataFrame) -> pd.DataFrame:
    station = _weather_station_for(result, exp_id)
    if not station:
        return pd.DataFrame()
    try:
        wth = resolve_dssat_paths(result.cfg)["weather"] / f"{station}.WTH"
        df = read_wth(wth)
    except Exception:
        return pd.DataFrame()
    if df.empty:
        return pd.DataFrame()
    if sim_long is not None and not sim_long.empty and "date" in sim_long.columns:
        dates = pd.to_datetime(sim_long["date"], errors="coerce").dropna()
        if not dates.empty:
            df = df[(df["date"] >= dates.min()) & (df["date"] <= dates.max())]
    value_cols = [c for c in ("TMIN", "TMAX", "SRAD") if c in df.columns]
    if not value_cols:
        return pd.DataFrame()
    out = df.melt(id_vars=["date"], value_vars=value_cols, var_name="variable", value_name="value")
    out["source_file"] = "WTH"
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    return out.dropna(subset=["value"])


def _treatments_for_plot(sim_long: pd.DataFrame, obs: pd.DataFrame, exp_id: str) -> list[int | None]:
    values = []
    if sim_long is not None and not sim_long.empty and "treatment" in sim_long.columns:
        values.extend(sim_long["treatment"].dropna().tolist())
    if obs is not None and not obs.empty and "treatment" in obs.columns:
        o = obs[obs["exp_id"].astype(str) == str(exp_id)]
        values.extend(o["treatment"].dropna().tolist())
    out = sorted({int(v) for v in values if pd.notna(v)})
    return out or [None]


def _filter_treatment(df: pd.DataFrame, treatment: int | None) -> pd.DataFrame:
    if treatment is None or df is None or df.empty or "treatment" not in df.columns:
        return df
    return df[pd.to_numeric(df["treatment"], errors="coerce") == int(treatment)].copy()


def _daily_series(g: pd.DataFrame, spec: dict) -> pd.DataFrame:
    if g.empty or not spec.get("daily_average"):
        return g
    if "date" not in g.columns or not g["date"].notna().any():
        return g
    group_cols = ["date"]
    if "variable" in g.columns:
        group_cols.append("variable")
    return g.groupby(group_cols, as_index=False)["value"].mean()


def _soil_water_nitrogen_frame(sim_long: pd.DataFrame) -> pd.DataFrame:
    if sim_long is None or sim_long.empty:
        return pd.DataFrame(columns=["date", "variable", "value"])
    frames = []
    source = sim_long.get("source_file", pd.Series("", index=sim_long.index)).astype(str)

    water = sim_long[
        source.isin(["SoilWat.OUT", "SoilWater.OUT"]) &
        sim_long["variable"].astype(str).isin(["SWTD", "SW"])
    ].copy()
    if not water.empty:
        water["variable"] = water["variable"].astype(str).map({"SWTD": "soil water", "SW": "soil water"})
        frames.append(water.groupby(["date", "variable"], as_index=False)["value"].mean())

    n = sim_long[
        source.isin(["SoilNi.OUT"]) &
        sim_long["variable"].astype(str).isin([
            "NIAD", "NITD", "NHTD", "NO3", "NH4", "SNO3", "SNH4",
            *[f"NT{i}D" for i in range(1, 21)],
        ])
    ].copy()
    if not n.empty:
        n["date"] = pd.to_datetime(n["date"], errors="coerce")
        n = n.dropna(subset=["date"])
        n = n.groupby(["date", "variable"], as_index=False)["value"].sum()
        wide = n.pivot_table(index="date", columns="variable", values="value", aggfunc="sum")
        if "NIAD" in wide.columns:
            total = wide["NIAD"].fillna(0.0)
            found = True
        else:
            total = pd.Series(0.0, index=wide.index)
            found = False
            for col in ("NITD", "NHTD", "NO3", "SNO3", "NH4", "SNH4"):
                if col in wide.columns:
                    total = total.add(wide[col].fillna(0.0), fill_value=0.0)
                    found = True
            if not found:
                layer_cols = [c for c in wide.columns if str(c).startswith("NT") and str(c).endswith("D")]
                for col in layer_cols:
                    total = total.add(wide[col].fillna(0.0), fill_value=0.0)
                    found = True
        if found:
            frames.append(pd.DataFrame({
                "date": total.index,
                "variable": "plant available N",
                "value": total.to_numpy(dtype=float),
            }))

    if not frames:
        return pd.DataFrame(columns=["date", "variable", "value"])
    return pd.concat(frames, ignore_index=True, sort=False).dropna(subset=["date", "value"])


def _plot_series_panel(ax, sim_long: pd.DataFrame, obs: pd.DataFrame, exp_id: str,
                       spec: dict, aliases: dict[str, set[str]], color_offset: int = 0,
                       treatment: int | None = None) -> bool:
    cmap = plt.get_cmap("tab10")
    plotted = False
    sources = set(spec.get("sources") or [])
    for i, (var, label) in enumerate(spec.get("series", [])):
        color = cmap((i + color_offset) % 10)
        g = sim_long[sim_long["variable"].astype(str) == var].copy()
        if sources and "source_file" in g.columns:
            g = g[g["source_file"].isin(sources)]
        if not g.empty:
            g = _daily_series(g, spec)
            xcol = "date" if "date" in g.columns and g["date"].notna().any() else "DAP"
            groups = g.groupby("treatment", dropna=False) if "treatment" in g.columns else [(None, g)]
            for _trt, gg in groups:
                gg = gg.sort_values(xcol)
                ax.plot(gg[xcol], gg["value"], color=color, lw=1.2, alpha=0.75, label=f"sim {label}")
                plotted = True
        o = _obs_for_variable(obs, exp_id, var, aliases, treatment=treatment)
        if not o.empty:
            ox = o["date"] if o["date"].notna().any() else pd.to_numeric(o.get("value"), errors="coerce")
            if o["date"].notna().any():
                ax.scatter(ox, o["value"], color=color, edgecolor="black", linewidth=0.35,
                           s=22, marker="o", label=f"obs {label}", zorder=5)
                plotted = True
    ax.set_title(spec["title"], fontsize=9, fontweight="bold")
    ax.tick_params(axis="x", labelrotation=30, labelsize=7)
    if plotted:
        handles, labels = ax.get_legend_handles_labels()
        dedup = {}
        for h, l in zip(handles, labels):
            dedup.setdefault(l, h)
        ax.legend(dedup.values(), dedup.keys(), fontsize=6, frameon=False, loc="best")
    else:
        ax.text(0.5, 0.5, "not available", ha="center", va="center",
                transform=ax.transAxes, color="#888780", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    return plotted


def _plot_soil_panel(ax, sim_long: pd.DataFrame, treatment: int | None = None) -> bool:
    cmap = plt.get_cmap("tab10")
    plotted = False
    soil = _soil_water_nitrogen_frame(sim_long)
    for i, label in enumerate(["soil water", "plant available N"]):
        g = soil[soil["variable"] == label].copy()
        if g.empty:
            continue
        g = g.sort_values("date")
        ax.plot(g["date"], g["value"], color=cmap(i), lw=1.2, alpha=0.8, label=f"sim {label}")
        plotted = True
    ax.set_title("Soil water and nitrogen", fontsize=9, fontweight="bold")
    ax.tick_params(axis="x", labelrotation=30, labelsize=7)
    if plotted:
        ax.legend(fontsize=6, frameon=False, loc="best")
    else:
        ax.text(0.5, 0.5, "not available", ha="center", va="center",
                transform=ax.transAxes, color="#888780", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    return plotted


def _planting_date_from_plantgro(pg: pd.DataFrame, treatment: int | None):
    if treatment is None or pg is None or pg.empty or "date" not in pg.columns or "DAP" not in pg.columns:
        return pd.NaT
    sub = _filter_treatment(pg, treatment)
    if sub.empty:
        return pd.NaT
    dates = pd.to_datetime(sub["date"], errors="coerce")
    dap = pd.to_numeric(sub["DAP"], errors="coerce")
    ok = dates.notna() & dap.notna()
    if not ok.any():
        return pd.NaT
    planted = dates[ok] - pd.to_timedelta(dap[ok], unit="D")
    return planted.iloc[0] if not planted.empty else pd.NaT


def _stage_observed_dap(rows: pd.DataFrame, pg: pd.DataFrame, treatment: int | None) -> pd.Series:
    if rows.empty:
        return pd.Series(dtype=float)
    planted = _planting_date_from_plantgro(pg, treatment)
    vals = []
    for _, row in rows.iterrows():
        date = pd.to_datetime(row.get("date"), errors="coerce")
        if pd.notna(date) and pd.notna(planted):
            vals.append(float((date.normalize() - planted.normalize()).days))
            continue
        raw = pd.to_numeric(row.get("value"), errors="coerce")
        if pd.notna(raw) and raw < 1000:
            vals.append(float(raw))
    return pd.Series(vals, dtype=float)


def _plot_stage_panel(ax, result, spawn, exp_id: str, treatment: int | None = None) -> bool:
    stage_defs = [
        ("EDAP", {"EDAT", "EDATE"}, "emergence"),
        ("ADAP", {"ADAT"}, "anthesis"),
        ("MDAP", {"MDAT"}, "maturity"),
        ("HDAP", {"HDAT"}, "harvest"),
    ]
    ylabels = []
    plotted = False
    ev = getattr(spawn, "evaluate", pd.DataFrame())
    obs = result.obs.table if hasattr(result, "obs") else pd.DataFrame()
    cmap = plt.get_cmap("tab10")
    pg = getattr(spawn, "plantgro", pd.DataFrame())
    for yi, (dap_var, date_vars, label) in enumerate(stage_defs):
        ylabels.append(label)
        if ev is not None and not ev.empty:
            g = ev[ev["variable"].astype(str) == dap_var]
            g = _filter_treatment(g, treatment)
            if not g.empty:
                vals = pd.to_numeric(g["sim"], errors="coerce").dropna()
                if not vals.empty:
                    ax.scatter(vals, [yi] * len(vals), marker="x", color=cmap(0), s=36,
                               label="simulated" if yi == 0 else None)
                    plotted = True
        if obs is not None and not obs.empty:
            o = obs[(obs["exp_id"].astype(str) == str(exp_id)) &
                    (obs["variable"].astype(str).isin(set(date_vars) | {dap_var, label}))].copy()
            o = _filter_treatment(o, treatment)
            if not o.empty:
                vals = _stage_observed_dap(o, pg, treatment).dropna()
                if not vals.empty:
                    ax.scatter(vals, [yi] * len(vals), marker="o", color=cmap(1),
                               edgecolor="black", linewidth=0.35, s=28,
                               label="observed" if yi == 0 else None)
                    plotted = True
    ax.set_title("Growth stages", fontsize=9, fontweight="bold")
    ax.set_xlabel("days after planting")
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=7)
    if plotted:
        handles, labels = ax.get_legend_handles_labels()
        labels_seen = {}
        for h, l in zip(handles, labels):
            if l:
                labels_seen.setdefault(l, h)
        ax.legend(labels_seen.values(), labels_seen.keys(), fontsize=6, frameon=False, loc="best")
    else:
        ax.text(0.5, 0.5, "not available", ha="center", va="center",
                transform=ax.transAxes, color="#888780", fontsize=9)
    return plotted


def plot_experiment_diagnostics(result, best_spawns, figdir) -> list[Path]:
    """Write one 3x3 diagnostic time-series panel per calibrated treatment.

    Panels cover LAI, stages, canopy dimensions, total/organ biomass, stress,
    soil water/nitrogen, weather, and tissue nitrogen when the corresponding
    DSSAT outputs/observations are available.
    """
    if not best_spawns:
        return []
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    obs = result.obs.table if hasattr(result, "obs") else pd.DataFrame()
    aliases = _obs_aliases(result)
    paths: list[Path] = []
    for exp_id, spawn in best_spawns.items():
        if getattr(spawn, "status", "") not in {"success", "cached"}:
            continue
        sim_long = _combined_output_long(spawn)
        for treatment in _treatments_for_plot(sim_long, obs, str(exp_id)):
            sim_trt = _filter_treatment(sim_long, treatment)
            weather = _weather_fallback_long(result, str(exp_id), sim_trt)
            if not weather.empty:
                sim_trt = pd.concat([sim_trt, weather], ignore_index=True, sort=False)
            fig, axes = plt.subplots(3, 3, figsize=(14, 10.5), squeeze=False)
            for idx, spec in enumerate(_EXPERIMENT_PANEL_SPECS):
                ax = axes[idx // 3][idx % 3]
                if spec.get("stage"):
                    _plot_stage_panel(ax, result, spawn, exp_id, treatment=treatment)
                elif spec.get("soil"):
                    _plot_soil_panel(ax, sim_trt, treatment=treatment)
                else:
                    _plot_series_panel(ax, sim_trt, obs, str(exp_id), spec, aliases,
                                       color_offset=idx, treatment=treatment)
            trt_label = f"T{int(treatment)}" if treatment is not None else "TNA"
            fig.suptitle(f"{exp_id} {trt_label}: observed vs simulated diagnostic time series",
                         fontsize=13, fontweight="bold")
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            path = figdir / f"fig_experiment_{exp_id}_{trt_label}_3x3.png"
            fig.savefig(path, dpi=140, bbox_inches="tight")
            plt.close(fig)
            paths.append(path)
    return paths


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


def plot_sensitivity(result, path):
    """Tornado bar of parameter influence from the screening stage (if it ran)."""
    ranking = (result.extras or {}).get("sensitivity")
    if ranking is None or ranking.empty:
        return
    metric = "mu_star" if "mu_star" in ranking.columns else "ST"
    r = ranking.sort_values(metric, ascending=True)   # smallest at bottom
    fig, ax = plt.subplots(figsize=(6.5, max(3, 0.4 * len(r) + 1)))
    ax.barh(r["parameter"], r[metric], color=_ACCENT)
    label = "mu*  (mean |elementary effect|)" if metric == "mu_star" else "ST  (total Sobol index)"
    ax.set_xlabel(label)
    ax.set_title("Parameter influence (higher = matters more)")
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_mcmc_trace(result, path):
    """MCMC only: log-posterior trace per walker (visual mixing/convergence check)."""
    chain = (result.extras or {}).get("mcmc_chain")
    if chain is None or chain.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 4))
    for _w, g in chain.groupby("walker"):
        ax.plot(g["step"], g["logpost"], lw=0.7, alpha=0.6)
    ax.set_xlabel("step"); ax.set_ylabel("log-posterior")
    acc = (result.extras or {}).get("acceptance")
    title = "MCMC log-posterior traces (walkers should overlap once mixed)"
    if acc is not None:
        title += f"  —  acceptance {acc:.2f}"
    ax.set_title(title)
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


def objective_breakdown_table(result) -> pd.DataFrame:
    """Return per-experiment/per-variable objective components for the best fit."""
    resid = getattr(result.best, "residuals", pd.DataFrame())
    if resid is None or resid.empty:
        return pd.DataFrame(columns=[
            "exp_id", "user_var", "kind", "n", "mean_loss", "weighted_loss",
            "RMSE", "MBE", "mean_obs", "mean_sim",
        ])
    df = resid.copy()
    if "_loss" not in df.columns:
        sigma_source = df["sigma"] if "sigma" in df.columns else pd.Series(1.0, index=df.index)
        sigma = pd.to_numeric(sigma_source, errors="coerce").replace(0, np.nan)
        df["_loss"] = (pd.to_numeric(df["resid"], errors="coerce") / sigma) ** 2
    if "weight" not in df.columns:
        df["weight"] = 1.0
    group_cols = [c for c in ("exp_id", "user_var", "kind") if c in df.columns]
    rows = []
    for keys, g in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        rec = dict(zip(group_cols, keys))
        resid_vals = pd.to_numeric(g["resid"], errors="coerce")
        rec.update({
            "n": int(len(g)),
            "mean_loss": float(pd.to_numeric(g["_loss"], errors="coerce").mean()),
            "weighted_loss": float((pd.to_numeric(g["_loss"], errors="coerce") *
                                    pd.to_numeric(g["weight"], errors="coerce")).sum()),
            "RMSE": float(np.sqrt(np.nanmean(resid_vals ** 2))),
            "MBE": float(resid_vals.mean()),
            "mean_obs": float(pd.to_numeric(g["obs"], errors="coerce").mean()),
            "mean_sim": float(pd.to_numeric(g["sim"], errors="coerce").mean()),
        })
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(group_cols).reset_index(drop=True)


def _date_str(value) -> str:
    ts = pd.to_datetime(value, errors="coerce")
    return "" if pd.isna(ts) else str(ts.date())


def _observed_stage_date(obs: pd.DataFrame, exp_id: str, treatment: int, variables: set[str]):
    if obs is None or obs.empty:
        return pd.NaT
    rows = obs[
        (obs["exp_id"].astype(str) == str(exp_id))
        & (pd.to_numeric(obs["treatment"], errors="coerce") == int(treatment))
        & (obs["variable"].astype(str).isin(variables))
    ]
    if rows.empty:
        return pd.NaT
    dates = pd.to_datetime(rows["date"], errors="coerce").dropna()
    if not dates.empty:
        return dates.iloc[0]
    return pd.NaT


def _header_token_spans(header: str) -> dict[str, tuple[int, int]]:
    spans = {}
    for match in re.finditer(r"\S+", header):
        spans[match.group().lstrip("@")] = (match.start(), match.end())
    return spans


def _selected_cultivar_name(cfg: dict, exp_id: str, treatment: int) -> str:
    try:
        crop = crop_for(cfg, (cfg.get("crops") or [{}])[0].get("code", "HM"))
        filex = Path(cfg["source"]["hemp_dir"]) / f"{exp_id}.{crop['filex_ext']}"
    except Exception:
        return ""
    if not filex.exists():
        return ""
    lines = filex.read_text(errors="replace").splitlines()

    cultivar_by_factor = {}
    try:
        start = next(i for i, ln in enumerate(lines) if ln.startswith("*CULTIVARS"))
    except StopIteration:
        start = None
    if start is not None:
        header = None
        for ln in lines[start + 1:]:
            if ln.startswith("*"):
                break
            if ln.lstrip().startswith("@"):
                header = ln.lstrip().lstrip("@").split()
                continue
            if header and ln.strip() and not ln.lstrip().startswith("!"):
                vals = ln.split(maxsplit=len(header) - 1)
                row = dict(zip(header, vals))
                factor = row.get("C")
                if factor:
                    name = row.get("CNAME") or row.get("INGENO") or ""
                    code = row.get("INGENO") or ""
                    cultivar_by_factor[str(int(float(factor)))] = str(name or code)

    try:
        start = next(i for i, ln in enumerate(lines) if ln.startswith("*TREATMENTS"))
    except StopIteration:
        return ""
    header = None
    spans = {}
    for ln in lines[start + 1:]:
        if ln.startswith("*"):
            break
        if ln.lstrip().startswith("@"):
            header = ln
            spans = _header_token_spans(ln)
            continue
        if not header or not ln.strip() or ln.lstrip().startswith("!"):
            continue
        try:
            trt = int(float(ln.split()[0]))
        except (ValueError, IndexError):
            continue
        if trt != int(treatment):
            continue
        if "CU" in spans:
            lo, hi = spans["CU"]
            try:
                factor = str(int(float(ln[lo:hi].strip())))
                return cultivar_by_factor.get(factor, factor)
            except ValueError:
                return ""
    return ""


def phenology_report_table(result, best_spawns=None) -> pd.DataFrame:
    """Return site-level planting, emergence, anthesis, and bias for best fit."""
    cols = [
        "site", "cultivar", "planting_date", "emergence_date", "observed_anthesis",
        "simulated_anthesis", "bias",
    ]
    resid = getattr(result.best, "residuals", pd.DataFrame())
    if resid is None or resid.empty:
        return pd.DataFrame(columns=cols)
    mask = pd.Series(False, index=resid.index)
    if "user_var" in resid.columns:
        mask = resid["user_var"].astype(str) == "anthesis"
    if "dssat" in resid.columns:
        mask = mask | (resid["dssat"].astype(str) == "ADAP")
    anth = resid[mask].copy()
    if anth.empty:
        return pd.DataFrame(columns=cols)

    obs = result.obs.table if hasattr(result, "obs") else pd.DataFrame()
    rows = []
    for _, row in anth.sort_values(["exp_id", "treatment"]).iterrows():
        exp_id = str(row["exp_id"])
        treatment = int(row["treatment"])
        spawn = (best_spawns or {}).get(exp_id)
        pg = getattr(spawn, "plantgro", pd.DataFrame()) if spawn is not None else pd.DataFrame()
        planted = _planting_date_from_plantgro(pg, treatment)
        emergence = _observed_stage_date(obs, exp_id, treatment, {"EDAT", "EDATE"})
        observed = _observed_stage_date(obs, exp_id, treatment, {"ADAT"})
        if pd.isna(observed) and pd.notna(planted):
            observed = planted + pd.to_timedelta(float(row["obs"]), unit="D")
        simulated = pd.NaT
        if pd.notna(planted) and pd.notna(row["sim"]):
            simulated = planted + pd.to_timedelta(float(row["sim"]), unit="D")
        bias = pd.to_numeric(row.get("resid"), errors="coerce")
        rows.append({
            "site": exp_id,
            "cultivar": _selected_cultivar_name(getattr(result, "cfg", {}) or {}, exp_id, treatment),
            "planting_date": _date_str(planted),
            "emergence_date": _date_str(emergence),
            "observed_anthesis": _date_str(observed),
            "simulated_anthesis": _date_str(simulated),
            "bias": int(round(float(bias))) if pd.notna(bias) else np.nan,
        })
    return pd.DataFrame(rows, columns=cols)


def sample_phenology_residuals_table(result) -> pd.DataFrame:
    """Return anthesis residuals for every evaluated sample/site.

    This is useful when the calibration target is not only the default scalar
    objective, but also balance across environments, such as minimizing the
    maximum absolute anthesis bias.
    """
    cols = [
        "sample_id", "site", "treatment", "user_var", "observed_dap",
        "simulated_dap", "bias", "abs_bias",
    ]
    rows = []
    for sample_id, ores in (getattr(result, "obj_results", {}) or {}).items():
        resid = getattr(ores, "residuals", pd.DataFrame())
        if resid is None or resid.empty:
            continue
        mask = pd.Series(False, index=resid.index)
        if "user_var" in resid.columns:
            mask = resid["user_var"].astype(str) == "anthesis"
        if "dssat" in resid.columns:
            mask = mask | (resid["dssat"].astype(str) == "ADAP")
        anth = resid[mask].copy()
        for _, row in anth.iterrows():
            bias = pd.to_numeric(row.get("resid"), errors="coerce")
            rows.append({
                "sample_id": sample_id,
                "site": str(row.get("exp_id", "")),
                "treatment": int(row["treatment"]) if pd.notna(row.get("treatment")) else np.nan,
                "user_var": str(row.get("user_var", "")),
                "observed_dap": float(row["obs"]) if pd.notna(row.get("obs")) else np.nan,
                "simulated_dap": float(row["sim"]) if pd.notna(row.get("sim")) else np.nan,
                "bias": float(bias) if pd.notna(bias) else np.nan,
                "abs_bias": abs(float(bias)) if pd.notna(bias) else np.nan,
            })
    return pd.DataFrame(rows, columns=cols)


def balanced_candidates_table(result) -> pd.DataFrame:
    """Rank samples by worst-site anthesis bias, then RMSE and mean bias."""
    phen = sample_phenology_residuals_table(result)
    if phen.empty:
        return pd.DataFrame(columns=[
            "sample_id", "n", "max_abs_bias", "RMSE", "MBE", "score",
        ])
    grouped = phen.groupby("sample_id", dropna=False)
    rows = []
    score_map = {}
    if hasattr(result, "design") and result.design is not None and not result.design.empty:
        if "sample_id" in result.design.columns and "score" in result.design.columns:
            score_map = dict(zip(result.design["sample_id"].astype(str), result.design["score"]))
    for sample_id, g in grouped:
        bias = pd.to_numeric(g["bias"], errors="coerce")
        rec = {
            "sample_id": sample_id,
            "n": int(bias.notna().sum()),
            "max_abs_bias": float(np.nanmax(np.abs(bias))) if bias.notna().any() else np.nan,
            "RMSE": float(np.sqrt(np.nanmean(bias ** 2))) if bias.notna().any() else np.nan,
            "MBE": float(np.nanmean(bias)) if bias.notna().any() else np.nan,
            "score": score_map.get(str(sample_id), np.nan),
        }
        for _, row in g.sort_values(["site", "treatment"]).iterrows():
            site = str(row["site"])
            rec[f"{site}_obs_dap"] = row["observed_dap"]
            rec[f"{site}_sim_dap"] = row["simulated_dap"]
            rec[f"{site}_bias"] = row["bias"]
        rows.append(rec)
    return (
        pd.DataFrame(rows)
        .sort_values(["max_abs_bias", "RMSE", "score"], na_position="last")
        .reset_index(drop=True)
    )


def spawn_manifest_table(result, best_spawns=None) -> pd.DataFrame:
    """Return durable spawn metadata from a calibration result."""
    manifest = (result.extras or {}).get("spawn_manifest")
    if manifest is not None and not manifest.empty:
        return manifest.copy()
    if not best_spawns:
        return pd.DataFrame(columns=[
            "sample_id", "exp_id", "theta_hash", "status", "message", "run_dir", "theta_json",
        ])
    theta = getattr(result, "best_theta", {}) or {}
    theta_jsonable = {k: (v.item() if hasattr(v, "item") else v) for k, v in theta.items()}
    rows = []
    for exp_id, res in best_spawns.items():
        rows.append({
            "sample_id": "best",
            "exp_id": exp_id,
            "theta_hash": theta_hash(theta) if theta else "",
            "status": getattr(res, "status", ""),
            "message": getattr(res, "message", ""),
            "run_dir": str(getattr(res, "run_dir", "")),
            "theta_json": json.dumps(theta_jsonable, sort_keys=True),
            **{f"theta_{k}": v for k, v in theta.items()},
        })
    return pd.DataFrame(rows)


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

    phenology = phenology_report_table(result, best_spawns=best_spawns)
    if not phenology.empty:
        phenology.to_csv(outdir / "phenology_report.csv", index=False)
        paths["phenology_report"] = outdir / "phenology_report.csv"

    objective_breakdown = objective_breakdown_table(result)
    objective_breakdown.round(6).to_csv(outdir / "objective_breakdown.csv", index=False)
    paths["objective_breakdown"] = outdir / "objective_breakdown.csv"

    sample_phenology = sample_phenology_residuals_table(result)
    if not sample_phenology.empty:
        sample_phenology.round(6).to_csv(outdir / "sample_phenology_residuals.csv", index=False)
        paths["sample_phenology_residuals"] = outdir / "sample_phenology_residuals.csv"
        balanced = balanced_candidates_table(result)
        balanced.round(6).to_csv(outdir / "balanced_candidates.csv", index=False)
        paths["balanced_candidates"] = outdir / "balanced_candidates.csv"

    spawn_manifest = spawn_manifest_table(result, best_spawns=best_spawns)
    if spawn_manifest is not None and not spawn_manifest.empty:
        spawn_manifest.to_csv(outdir / "manifest.csv", index=False)
        (outdir / "manifest.json").write_text(
            json.dumps(spawn_manifest.to_dict(orient="records"), indent=2, default=str),
            encoding="utf-8",
        )
        paths["manifest"] = outdir / "manifest.csv"
        paths["manifest_json"] = outdir / "manifest.json"

    (outdir / "best_theta.json").write_text(json.dumps(result.best_theta, indent=2))
    if result.glue is not None:
        from .engines.glue import posterior_summary
        posterior_summary(result.glue, result.space.names).round(4).to_csv(
            outdir / "posterior_summary.csv", index=False)
        paths["posterior_summary"] = outdir / "posterior_summary.csv"

    if not result.best.per_exp_var.empty:
        result.best.per_exp_var.round(3).to_csv(outdir / "fit_by_experiment.csv", index=False)

    # Extra data tables for the optional stages, when they ran.
    sens = (result.extras or {}).get("sensitivity")
    if sens is not None and not sens.empty:
        sens.round(4).to_csv(outdir / "sensitivity_ranking.csv", index=False)
        paths["sensitivity_ranking"] = outdir / "sensitivity_ranking.csv"
    chain = (result.extras or {}).get("mcmc_chain")
    if chain is not None and not chain.empty:
        chain.round(4).to_csv(outdir / "mcmc_chain.csv", index=False)
        paths["mcmc_chain"] = outdir / "mcmc_chain.csv"
    optimizer_history = (result.extras or {}).get("optimizer_history")
    if optimizer_history:
        history = pd.DataFrame(optimizer_history)
        history.to_csv(outdir / "optimizer_history.csv", index=False)
        paths["optimizer_history"] = outdir / "optimizer_history.csv"

    plot_param_posteriors(result, figdir / "fig_param_posteriors.png")
    plot_score_funnel(result, figdir / "fig_score_funnel.png")
    plot_ess_trajectory(result, figdir / "fig_ess_trajectory.png")
    plot_mcmc_trace(result, figdir / "fig_mcmc_trace.png")
    plot_sensitivity(result, figdir / "fig_sensitivity.png")
    plot_obs_vs_sim(result, figdir / "fig_obs_vs_sim.png")
    plot_obs_vs_sim_by_category(result, figdir / "fig_obs_vs_sim_by_category.png")
    plot_fit_bars(result, figdir / "fig_fit_bars.png")
    if best_spawns:
        plot_timeseries(result, best_spawns, figdir / "fig_timeseries.png")
        experiment_panels = plot_experiment_diagnostics(result, best_spawns, figdir)
        if experiment_panels:
            paths["experiment_diagnostics"] = experiment_panels

    if result.nsga2 is not None:
        _plot_pareto(result.nsga2, figdir / "fig_pareto.png")
        result.nsga2.front().round(4).to_csv(outdir / "pareto_front.csv", index=False)

    for key in ("param_posteriors", "score_funnel", "ess_trajectory", "mcmc_trace",
                "sensitivity", "obs_vs_sim", "obs_vs_sim_by_category", "fit_bars",
                "timeseries", "pareto"):
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
