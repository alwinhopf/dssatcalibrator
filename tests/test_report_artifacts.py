from types import SimpleNamespace

import pandas as pd

from dssatcalibrator import viz


def test_make_report_writes_manifest_and_objective_breakdown(tmp_path, monkeypatch):
    for name in (
        "plot_param_posteriors",
        "plot_score_funnel",
        "plot_ess_trajectory",
        "plot_mcmc_trace",
        "plot_sensitivity",
        "plot_obs_vs_sim",
        "plot_obs_vs_sim_by_category",
        "plot_fit_bars",
    ):
        monkeypatch.setattr(viz, name, lambda *args, **kwargs: None)

    residuals = pd.DataFrame({
        "exp_id": ["E1", "E1", "E2"],
        "user_var": ["biomass", "biomass", "yield"],
        "kind": ["timeseries", "timeseries", "scalar"],
        "obs": [10.0, 12.0, 100.0],
        "sim": [11.0, 13.0, 90.0],
        "resid": [1.0, 1.0, -10.0],
        "sigma": [1.0, 1.0, 10.0],
        "weight": [1.0, 1.0, 2.0],
        "_loss": [1.0, 1.0, 1.0],
    })
    best = SimpleNamespace(
        residuals=residuals,
        per_var={
            "biomass": {"n": 2, "RMSE": 1.0, "nRMSE_pct": 9.1, "MBE": 1.0, "d": 0.9, "EF": 0.5, "R2": 1.0},
            "yield": {"n": 1, "RMSE": 10.0, "nRMSE_pct": 10.0, "MBE": -10.0, "d": 0.8, "EF": 0.0, "R2": float("nan")},
        },
        per_exp_var=pd.DataFrame(),
    )
    result = SimpleNamespace(
        cfg={"method": {"bayesian": {"behavioural_quantile": 0.1}}},
        space=SimpleNamespace(names=["P1"], low=[0.0], high=[2.0], start=[1.0]),
        design=pd.DataFrame({"sample_id": [0], "P1": [1.0], "score": [1.0], "loglik": [-1.0], "n_obs": [3], "weight": [1.0]}),
        best_theta={"P1": 1.0},
        best=best,
        glue=None,
        nsga2=None,
        extras={
            "spawn_manifest": pd.DataFrame({
                "sample_id": [0, 0],
                "exp_id": ["E1", "E2"],
                "status": ["success", "success"],
                "run_dir": ["runs/E1", "runs/E2"],
                "theta_json": ['{"P1": 1.0}', '{"P1": 1.0}'],
            })
        },
    )

    paths = viz.make_report(result, tmp_path)

    assert paths["manifest"].name == "manifest.csv"
    assert (tmp_path / "manifest.json").exists()
    assert (tmp_path / "objective_breakdown.csv").exists()
    manifest = pd.read_csv(tmp_path / "manifest.csv")
    breakdown = pd.read_csv(tmp_path / "objective_breakdown.csv")
    assert manifest.shape[0] == 2
    assert {"exp_id", "user_var", "weighted_loss", "RMSE"} <= set(breakdown.columns)
    assert breakdown.query("exp_id == 'E2'")["weighted_loss"].iloc[0] == 2.0
