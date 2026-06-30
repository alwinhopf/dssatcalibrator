from dssatcalibrator.config import active_parameters, validate_config
from dssatcalibrator.spaces import ParameterSpace, expand_parameter_specs
from dssatcalibrator.spawn import _effective_theta, _partition_theta


def pooled_cfg():
    return {
        "experiments": ["E1", "E2"],
        "parameters": {
            "genetic_species": {
                "SPEP": {"active": True, "min": 1.0, "max": 3.0, "start": 2.0, "spe_key": "SPEP"},
            },
            "genetic_cultivar": {
                "CULP": {"active": True, "min": 10.0, "max": 20.0, "start": 15.0, "scope": "experiment"},
            },
            "genetic_ecotype": {
                "ECOP": {"active": True, "min": 100.0, "max": 200.0, "start": 150.0, "pooling": "per_experiment"},
            },
        },
    }


def test_experiment_scoped_parameters_expand_per_experiment():
    cfg = pooled_cfg()
    specs = expand_parameter_specs(cfg, active_parameters(cfg))

    assert [s["name"] for s in specs] == [
        "SPEP",
        "CULP__E1",
        "CULP__E2",
        "ECOP__E1",
        "ECOP__E2",
    ]
    assert [s["base_name"] for s in specs] == ["SPEP", "CULP", "CULP", "ECOP", "ECOP"]
    assert [s["scope"] for s in specs] == ["global", "experiment", "experiment", "experiment", "experiment"]
    assert [s.get("exp_id") for s in specs] == [None, "E1", "E2", "E1", "E2"]

    space = ParameterSpace.from_config(cfg)
    assert space.names == ["SPEP", "CULP__E1", "CULP__E2", "ECOP__E1", "ECOP__E2"]
    assert space.start.tolist() == [2.0, 15.0, 15.0, 150.0, 150.0]


def test_spawn_partition_uses_only_values_for_current_experiment():
    cfg = pooled_cfg()
    specs = expand_parameter_specs(cfg, active_parameters(cfg))
    theta = {
        "SPEP": 2.5,
        "CULP__E1": 11.0,
        "CULP__E2": 19.0,
        "ECOP__E1": 111.0,
        "ECOP__E2": 199.0,
    }

    groups = _partition_theta(theta, specs, exp_id="E1")
    assert groups["genetic_species"] == {"SPEP": 2.5}
    assert groups["genetic_cultivar"] == {"CULP": 11.0}
    assert groups["genetic_ecotype"] == {"ECOP": 111.0}

    effective = _effective_theta(theta, specs, exp_id="E1")
    assert effective == {"SPEP": 2.5, "CULP": 11.0, "ECOP": 111.0}

    e2_groups = _partition_theta(theta, specs, exp_id="E2")
    assert e2_groups["genetic_cultivar"] == {"CULP": 19.0}
    assert e2_groups["genetic_ecotype"] == {"ECOP": 199.0}


def test_validate_config_accepts_parameter_scope_and_rejects_unknown_scope():
    cfg = pooled_cfg()
    assert validate_config(cfg) is cfg

    bad = pooled_cfg()
    bad["parameters"]["genetic_cultivar"]["CULP"]["scope"] = "plot"
    try:
        validate_config(bad)
    except ValueError as exc:
        assert "scope 'plot'" in str(exc)
    else:
        raise AssertionError("validate_config accepted an unknown parameter scope")
