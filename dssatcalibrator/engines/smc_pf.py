"""SMC Particle Filter engine with Metropolis-Hastings mutation (preset A).

Sequentially assimilates time-ordered observations (time-series, followed by
end-of-season scalars/phenology), calculates particle weights based on the cumulative
likelihood, and triggers systematic resampling and parameter mutation (MCMC) when the
effective sample size (ESS) drops below the threshold.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
import pandas as pd

from .. import objective as obj
from ..config import resolve_exe
from ..runner import resolve_cores, run_many
from ..samplers import sample
from ..spaces import ParameterSpace


@dataclass
class SmcResult:
    design: pd.DataFrame          # design + a "weight" column (posterior weights)
    behavioural: pd.DataFrame     # the behavioural subset
    best_theta: dict
    best_sample_id: int
    threshold: float
    ess: float                    # effective sample size of the posterior weights
    obj_results: dict             # sid -> ObjectiveResult
    best: obj.ObjectiveResult
    ess_trace: list = field(default_factory=list)   # [{step,label,ess,n,resampled}, ...]
    initial_design: pd.DataFrame = None             # the initial ensemble (the prior)


def run_smc_pf(cfg: dict, progress: bool = True) -> SmcResult:
    # 1. Setup
    from ..orchestrator import _setup
    space, crop, exe, specs, run_root, obs, experiments, treatments = _setup(cfg)
    n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))

    method = cfg.get("method", {})
    bayesian_cfg = method.get("bayesian", {})
    
    # Read SMC parameters from config
    n_particles = int(bayesian_cfg.get("n_particles", 200))
    ess_frac = float(bayesian_cfg.get("ess_frac", 0.5))
    mutation_scale = float(bayesian_cfg.get("mutation_scale", 0.02))   # fixed-kernel sd (frac of range)
    move_kernel = str(bayesian_cfg.get("move_kernel", "adaptive"))     # "adaptive" | "fixed"
    kernel_floor = float(bayesian_cfg.get("kernel_floor", 0.01))       # min jitter (frac of range)
    kernel_scale = bayesian_cfg.get("kernel_scale", None)             # None -> optimal 2.38/sqrt(d)
    seed = int(cfg["calibrator"].get("seed", 42))

    # Dedicated RNG (do not pollute the global np.random state used elsewhere)
    rng = np.random.default_rng(seed)

    # 2. Sample initial design
    sample_engine = method.get("sample", {}).get("engine", "lhs")
    samples = sample(space, n=n_particles, engine=sample_engine, seed=seed, include_start=True)
    # Update n_particles if it got adjusted by sampler (some samplers round n)
    n_particles = len(samples)

    # Keep the initial ensemble as the "prior" for prior-vs-posterior plots.
    initial_design = pd.DataFrame(
        [{"sample_id": sid, **space.to_theta(samples.loc[sid].to_numpy())}
         for sid in range(n_particles)]
    )

    # 3. Run initial ensemble
    if progress:
        print(f"Running initial ensemble of {n_particles} particles...", flush=True)

    jobs = []
    idx = []
    for sid, row in samples.iterrows():
        theta = space.to_theta(row.to_numpy())
        for exp in experiments:
            jobs.append(dict(theta=dict(theta), exp_id=exp, cfg=cfg, crop=crop, param_specs=specs,
                             run_root=run_root, treatments=treatments[exp], exe=exe))
            idx.append((sid, exp))

    done = {"n": 0}
    total_jobs = len(jobs)
    def _cb(_res):
        done["n"] += 1
        if progress and (done["n"] % max(1, total_jobs // 20) == 0 or done["n"] == total_jobs):
            print(f"  initial spawns {done['n']}/{total_jobs}", flush=True)

    results = run_many(jobs, n_workers=n_workers, on_done=_cb if progress else None)

    # Group initial results by particle
    particles = []
    for sid in range(n_particles):
        p_results = {}
        for (job_sid, exp), res in zip(idx, results):
            if job_sid == sid:
                p_results[exp] = res
        
        theta = space.to_theta(samples.loc[sid].to_numpy())
        resid_df = obj.build_residuals(p_results, obs.table, cfg)
        
        particles.append({
            "theta": theta,
            "results": p_results,
            "residuals": resid_df,
            "loglik": -1e10 if resid_df.empty else 0.0,
        })

    # 4. Define the sequential assimilation steps from the matched observations.
    # A step is a dict {"date", "var", "label"}:
    #   * time-series  -> one step per unique date  (var=None: all variables on that date)
    #   * scalars/phen -> one step PER variable     (date=NaT, var=user_var)
    # Splitting the end-of-season scalars per variable stops a single huge "scalars"
    # step from dominating the filter and collapsing ESS in one shot.
    unique_dates, scalar_vars = set(), set()
    for p in particles:
        if p["residuals"].empty:
            continue
        for _, r in p["residuals"].iterrows():
            if pd.isna(r["date"]):
                scalar_vars.add(r["user_var"])
            else:
                unique_dates.add(r["date"])

    steps = [{"date": d, "var": None, "label": pd.Timestamp(d).strftime("%Y-%m-%d")}
             for d in sorted(unique_dates)]
    steps += [{"date": pd.NaT, "var": v, "label": f"scalar:{v}"}
              for v in sorted(scalar_vars)]

    if progress:
        print(f"Identified {len(steps)} sequential assimilation steps "
              f"({len(unique_dates)} time-series dates, {len(scalar_vars)} end-of-season "
              f"scalar variables).", flush=True)

    # Helper: log-likelihood contribution of one step for a residual table.
    def compute_loglik_contrib(resid_df: pd.DataFrame, step: dict) -> float:
        if resid_df.empty:
            return -1e10
        if step["var"] is None:                      # a time-series date
            mask = resid_df["date"] == step["date"]
        else:                                        # a single end-of-season scalar
            mask = resid_df["date"].isna() & (resid_df["user_var"] == step["var"])
        sub = resid_df[mask]
        if sub.empty:
            return 0.0
        return float(-0.5 * np.sum(((sub["resid"] / sub["sigma"]) ** 2) * sub["weight"]))

    # Accumulated log-likelihoods for each particle
    loglik_accum = np.array([p["loglik"] for p in particles], dtype=float)

    n_dim = len(space.names)
    default_c = 2.38 / np.sqrt(max(n_dim, 1))        # optimal random-walk scaling

    def move_sd_by_name() -> dict[str, float]:
        """Per-parameter proposal sd for the MH move.

        Adaptive: scale to the current weighted ensemble spread (with a floor so the
        kernel never collapses to zero -> fights particle impoverishment).
        Fixed: a constant fraction of each parameter's range.
        """
        sd = {}
        if move_kernel == "adaptive":
            mat = np.array([[p["theta"][name] for name in space.names] for p in particles])
            wmean = np.average(mat, axis=0, weights=weights)
            wstd = np.sqrt(np.average((mat - wmean) ** 2, axis=0, weights=weights))
            c = float(kernel_scale) if kernel_scale is not None else default_c
            for j, name in enumerate(space.names):
                rng_j = space.high[j] - space.low[j]
                sd[name] = float(max(c * wstd[j], kernel_floor * rng_j))
        else:
            for j, name in enumerate(space.names):
                sd[name] = float(mutation_scale * (space.high[j] - space.low[j]))
        return sd

    def perturb_theta(parent_theta: dict[str, float], sd: dict[str, float]) -> dict[str, float]:
        mutated = {}
        for name, val in parent_theta.items():
            spec = next(s for s in space.specs if s["name"] == name)
            low, high = float(spec["min"]), float(spec["max"])
            mutated[name] = float(np.clip(val + rng.normal(0.0, sd[name]), low, high))
        return mutated

    # 5. Sequential Assimilation Loop
    # Defaults so the post-process is well-defined even if there are no steps.
    weights = np.ones(n_particles) / n_particles
    ess_trace: list = []
    for step_idx, step in enumerate(steps):
        # 5.1. Update log-likelihoods with new observations
        for i in range(n_particles):
            if loglik_accum[i] > -1e9:
                contrib = compute_loglik_contrib(particles[i]["residuals"], step)
                loglik_accum[i] += contrib

        # 5.2. Compute normalized weights
        max_loglik = np.max(loglik_accum)
        if max_loglik == -1e10 or not np.isfinite(max_loglik):
            weights = np.ones(n_particles) / n_particles
        else:
            w = np.exp(loglik_accum - max_loglik)
            w[~np.isfinite(w)] = 0.0
            total = w.sum()
            weights = w / total if total > 0 else np.ones(n_particles) / n_particles

        # 5.3. Compute ESS
        ess = 1.0 / np.sum(weights ** 2)
        step_label = step["label"]
        ess_rec = {"step": step_idx + 1, "label": step_label, "ess": float(ess),
                   "n": n_particles, "resampled": False}
        ess_trace.append(ess_rec)
        if progress:
            print(f"Step {step_idx+1}/{len(steps)} ({step_label}) -> ESS: {ess:.1f}/{n_particles}", flush=True)

        # 5.4. Check Resampling Condition (do not resample on the last step)
        if ess < n_particles * ess_frac and step_idx < len(steps) - 1:
            ess_rec["resampled"] = True
            if progress:
                print(f"  ESS is below threshold ({n_particles * ess_frac:.1f}). Resampling and mutating...", flush=True)

            # Proposal sd from the CURRENT ensemble spread (before resampling)
            sd = move_sd_by_name()
            if progress:
                print(f"  move kernel ({move_kernel}) sd: "
                      + ", ".join(f"{k}={v:.3g}" for k, v in sd.items()), flush=True)

            # Resample particle indices with replacement using weights
            resampled_idx = rng.choice(n_particles, size=n_particles, p=weights)

            # Generate mutated parameter sets
            mutated_thetas = []
            for i in range(n_particles):
                parent_idx = resampled_idx[i]
                parent_theta = particles[parent_idx]["theta"]
                mutated_thetas.append(perturb_theta(parent_theta, sd))

            # Run simulations for mutated particles in parallel
            mutation_jobs = []
            mut_idx = []
            for i in range(n_particles):
                mutated_theta = mutated_thetas[i]
                for exp in experiments:
                    mutation_jobs.append(dict(
                        theta=dict(mutated_theta), exp_id=exp, cfg=cfg, crop=crop, param_specs=specs,
                        run_root=run_root, treatments=treatments[exp], exe=exe
                    ))
                    mut_idx.append((i, exp))

            mut_done = {"n": 0}
            mut_total = len(mutation_jobs)
            def _mut_cb(_res):
                mut_done["n"] += 1
                if progress and (mut_done["n"] % max(1, mut_total // 20) == 0 or mut_done["n"] == mut_total):
                    print(f"    mutation spawns {mut_done['n']}/{mut_total}", flush=True)

            mut_results = run_many(mutation_jobs, n_workers=n_workers, on_done=_mut_cb if progress else None)

            # Group mutation results and perform Metropolis-Hastings acceptance check
            new_particles = []
            accepted_count = 0

            for i in range(n_particles):
                parent_idx = resampled_idx[i]
                parent_particle = particles[parent_idx]

                # Extract simulation results for mutated particle i
                p_mut_results = {}
                for (job_i, exp), res in zip(mut_idx, mut_results):
                    if job_i == i:
                        p_mut_results[exp] = res
                
                # Check if all runs succeeded
                succeeded = all(res.status in ("success", "cached") for res in p_mut_results.values())
                
                accepted = False
                if succeeded:
                    mut_resids = obj.build_residuals(p_mut_results, obs.table, cfg)
                    # Move targets the posterior given ALL data up to the current step,
                    # so both candidate and parent must be scored over the *full* history
                    # 0..step_idx. (loglik_accum is reset at each resample and only holds
                    # the incremental likelihood since the last move, so it cannot be used
                    # here — doing so biases the MH ratio after the second resample.)
                    loglik_mutated = 0.0
                    loglik_parent = 0.0
                    for s_idx in range(step_idx + 1):
                        loglik_mutated += compute_loglik_contrib(mut_resids, steps[s_idx])
                        loglik_parent += compute_loglik_contrib(parent_particle["residuals"], steps[s_idx])

                    # Metropolis-Hastings acceptance ratio
                    alpha = np.exp(loglik_mutated - loglik_parent)
                    if not np.isnan(alpha) and rng.uniform(0, 1) < alpha:
                        accepted = True
                        new_particles.append({
                            "theta": mutated_thetas[i],
                            "results": p_mut_results,
                            "residuals": mut_resids,
                            "loglik": loglik_mutated,
                        })
                        accepted_count += 1
                
                if not accepted:
                    # Keep parent particle
                    new_particles.append({
                        "theta": parent_particle["theta"].copy(),
                        "results": parent_particle["results"],
                        "residuals": parent_particle["residuals"].copy(),
                        "loglik": loglik_accum[parent_idx],
                    })

            particles = new_particles
            # Reset weights to uniform since we resampled and mutated
            loglik_accum = np.zeros(n_particles, dtype=float)
            if progress:
                print(f"  Mutation complete. Accepted: {accepted_count}/{n_particles} ({100*accepted_count/n_particles:.1f}%)", flush=True)

    # 6. Post-process final ensemble: compute final total scores and logliks
    obj_results = {}
    rows = []
    
    if progress:
        print("Sequential assimilation complete. Computing final total scores...", flush=True)

    for i in range(n_particles):
        p = particles[i]
        o = obj.score(p["results"], obs.table, cfg)
        obj_results[i] = o
        rec = {"sample_id": i, **p["theta"],
               "score": o.score, "loglik": o.loglik, "n_obs": len(o.residuals)}
        rows.append(rec)

    design = pd.DataFrame(rows)
    design["weight"] = weights

    best_sample_id = int(design["score"].idxmin()) if not design.empty else 0
    best_theta = particles[best_sample_id]["theta"]
    best = obj_results[best_sample_id]

    q = float(bayesian_cfg.get("behavioural_quantile", 0.1))
    valid = design[np.isfinite(design["score"])]
    threshold = float(valid["score"].quantile(q)) if not valid.empty else float("inf")
    behavioural = design[design["score"] <= threshold].copy()

    final_ess = 1.0 / np.sum(weights ** 2)

    return SmcResult(
        design=design,
        behavioural=behavioural,
        best_theta=best_theta,
        best_sample_id=best_sample_id,
        threshold=threshold,
        ess=final_ess,
        obj_results=obj_results,
        best=best,
        ess_trace=ess_trace,
        initial_design=initial_design,
    )
