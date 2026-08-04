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
        "plot_timeseries",
        "plot_experiment_diagnostics",
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
            }),
            "optimizer_history": [
                {"iter": 16, "score": 2.0},
                {"iter": 32, "score": 1.0},
            ],
        },
    )

    paths = viz.make_report(result, tmp_path)

    assert paths["manifest"].name == "manifest.csv"
    assert (tmp_path / "manifest.json").exists()
    assert (tmp_path / "objective_breakdown.csv").exists()
    assert paths["optimizer_history"].name == "optimizer_history.csv"
    history = pd.read_csv(paths["optimizer_history"])
    assert history.to_dict(orient="records") == [
        {"iter": 16, "score": 2.0},
        {"iter": 32, "score": 1.0},
    ]
    manifest = pd.read_csv(tmp_path / "manifest.csv")
    breakdown = pd.read_csv(tmp_path / "objective_breakdown.csv")
    assert manifest.shape[0] == 2
    assert {"exp_id", "user_var", "weighted_loss", "RMSE"} <= set(breakdown.columns)
    assert breakdown.query("exp_id == 'E2'")["weighted_loss"].iloc[0] == 2.0


def test_make_report_falls_back_to_best_spawn_manifest(tmp_path, monkeypatch):
    for name in (
        "plot_param_posteriors",
        "plot_score_funnel",
        "plot_ess_trajectory",
        "plot_mcmc_trace",
        "plot_sensitivity",
        "plot_obs_vs_sim",
        "plot_obs_vs_sim_by_category",
        "plot_fit_bars",
        "plot_timeseries",
        "plot_experiment_diagnostics",
    ):
        monkeypatch.setattr(viz, name, lambda *args, **kwargs: None)

    best = SimpleNamespace(
        residuals=pd.DataFrame(),
        per_var={},
        per_exp_var=pd.DataFrame(),
    )
    result = SimpleNamespace(
        cfg={"method": {"bayesian": {"behavioural_quantile": 0.1}}},
        space=SimpleNamespace(names=["P1"], low=[0.0], high=[2.0], start=[1.0]),
        design=pd.DataFrame({"sample_id": [0], "P1": [1.25], "score": [1.0], "loglik": [-1.0], "n_obs": [1]}),
        best_theta={"P1": 1.25},
        best=best,
        glue=None,
        nsga2=None,
        extras={},
    )
    best_spawns = {
        "E1": SimpleNamespace(status="success", message="", run_dir="runs/E1"),
        "E2": SimpleNamespace(status="success", message="", run_dir="runs/E2"),
    }

    paths = viz.make_report(result, tmp_path, best_spawns=best_spawns)

    manifest = pd.read_csv(paths["manifest"])
    assert manifest["sample_id"].tolist() == ["best", "best"]
    assert manifest["exp_id"].tolist() == ["E1", "E2"]
    assert manifest["status"].tolist() == ["success", "success"]
    assert {
        "theta_hash", "full_theta_hash", "effective_theta_hash",
        "theta_json", "effective_theta_json", "theta_P1",
    } <= set(manifest.columns)
    assert (tmp_path / "manifest.json").exists()


def test_phenology_report_table_dates_and_bias(tmp_path, monkeypatch):
    for name in (
        "plot_param_posteriors",
        "plot_score_funnel",
        "plot_ess_trajectory",
        "plot_mcmc_trace",
        "plot_sensitivity",
        "plot_obs_vs_sim",
        "plot_obs_vs_sim_by_category",
        "plot_fit_bars",
        "plot_timeseries",
        "plot_experiment_diagnostics",
    ):
        monkeypatch.setattr(viz, name, lambda *args, **kwargs: None)

    best = SimpleNamespace(
        residuals=pd.DataFrame({
            "exp_id": ["E1"],
            "treatment": [1],
            "user_var": ["anthesis"],
            "dssat": ["ADAP"],
            "kind": ["phenology"],
            "obs": [40.0],
            "sim": [43.0],
            "resid": [3.0],
        }),
        per_var={},
        per_exp_var=pd.DataFrame(),
    )
    obs = pd.DataFrame({
        "exp_id": ["E1", "E1"],
        "treatment": [1, 1],
        "variable": ["EDAT", "ADAT"],
        "kind": ["phenology", "phenology"],
        "date": pd.to_datetime(["2021-04-10", "2021-05-10"]),
        "value": [21100.0, 21130.0],
        "sigma": [float("nan"), float("nan")],
        "weight": [1.0, 1.0],
    })
    result = SimpleNamespace(
        cfg={},
        space=SimpleNamespace(names=[]),
        design=pd.DataFrame({"sample_id": [0]}),
        best_theta={},
        best=best,
        glue=None,
        nsga2=None,
        extras={},
        obs=SimpleNamespace(table=obs),
    )
    plantgro = pd.DataFrame({
        "date": pd.to_datetime(["2021-04-01", "2021-04-02"]),
        "DAP": [0, 1],
        "treatment": [1, 1],
    })
    best_spawns = {"E1": SimpleNamespace(plantgro=plantgro)}

    paths = viz.make_report(result, tmp_path, best_spawns=best_spawns)

    table = pd.read_csv(paths["phenology_report"])
    assert table.columns.tolist() == [
        "site", "cultivar", "planting_date", "emergence_date", "observed_anthesis",
        "simulated_anthesis", "bias",
    ]
    direct = viz.phenology_report_table(result, best_spawns=best_spawns)
    assert direct.iloc[0].to_dict() == {
        "site": "E1",
        "cultivar": "",
        "planting_date": "2021-04-01",
        "emergence_date": "2021-04-10",
        "observed_anthesis": "2021-05-10",
        "simulated_anthesis": "2021-05-14",
        "bias": 3,
    }


def test_plot_experiment_diagnostics_writes_panel(tmp_path):
    dates = pd.date_range("2021-06-01", periods=3)
    plantgro = pd.DataFrame({
        "date": dates,
        "treatment": [1, 1, 1],
        "DAP": [1, 2, 3],
        "LAID": [0.2, 0.5, 0.8],
        "CHTD": [0.1, 0.2, 0.3],
        "CWAD": [10.0, 25.0, 40.0],
        "SWAD": [4.0, 10.0, 16.0],
        "LWAD": [3.0, 8.0, 12.0],
        "RWAD": [1.0, 2.0, 3.0],
        "GWAD": [0.0, 0.0, 1.0],
        "WSPD": [1.0, 0.9, 0.8],
        "NSTD": [1.0, 1.0, 0.9],
    })
    extra_long = pd.DataFrame({
        "source_file": ["Weather.OUT", "Weather.OUT", "Weather.OUT", "SoilWat.OUT", "SoilNi.OUT", "PlantN.OUT"],
        "treatment": [1, 1, 1, 1, 1, 1],
        "date": [dates[0], dates[0], dates[0], dates[1], dates[1], dates[2]],
        "variable": ["TMIN", "TMAX", "SRAD", "SWTD", "NIAD", "LN%D"],
        "value": [20.0, 30.0, 18.5, 120.0, 10.0, 3.2],
    })
    obs = pd.DataFrame({
        "exp_id": ["E1", "E1", "E1"],
        "treatment": [1, 1, 1],
        "variable": ["LAID", "CWAD", "LN%D"],
        "kind": ["timeseries", "timeseries", "timeseries"],
        "date": [dates[1], dates[2], dates[2]],
        "value": [0.45, 42.0, 3.1],
        "sigma": [float("nan")] * 3,
        "weight": [1.0] * 3,
    })
    result = SimpleNamespace(
        cfg={"engine": {"timeseries_outputs": {"LAI": "LAID", "biomass": "CWAD"}}},
        obs=SimpleNamespace(table=obs),
    )
    spawn = SimpleNamespace(
        status="success",
        plantgro=plantgro,
        evaluate=pd.DataFrame({
            "treatment": [1],
            "run": [1],
            "variable": ["ADAP"],
            "sim": [42.0],
            "meas": [40.0],
        }),
        outputs={"long": extra_long},
    )

    paths = viz.plot_experiment_diagnostics(result, {"E1": spawn}, tmp_path)

    assert len(paths) == 1
    assert paths[0].name == "fig_experiment_E1_T1_3x3.png"
    assert paths[0].exists()
    assert paths[0].stat().st_size > 0
