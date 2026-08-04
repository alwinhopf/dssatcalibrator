"""Tests for config loading and active-parameter enumeration."""
from copy import deepcopy
from pathlib import Path

import pytest

from dssatcalibrator import config as cfgmod
from dssatcalibrator.runner import resolve_cores
from dssatcalibrator.samplers import sample
from dssatcalibrator.spaces import ParameterSpace, expand_parameter_specs

REPO = Path(__file__).resolve().parents[1]
HEMP_CFG = REPO / "config_hemp.yaml"


def test_validate_config_accepts_reference_config():
    # The shipped reference config must always pass validation.
    cfg = cfgmod.load_config(HEMP_CFG)            # load_config validates by default
    assert cfgmod.validate_config(cfg) is cfg


def test_validate_config_accepts_custom_pipeline():
    cfg = deepcopy(cfgmod.load_config(HEMP_CFG))
    cfg["method"]["preset"] = "custom"
    assert cfgmod.validate_config(cfg) is cfg


def test_validate_config_collects_all_problems():
    bad = {
        "method": {"preset": "Z", "validation": {"scheme": "weekly"}},
        "objective": {"weighting": "bogus"},
        "execution": {"backend": "spark"},
        "assimilation": {"mode": "kalman"},
        "gating": {"cultivar": "loose"},
        "calibrator": {"num_cores": -2},
        "parameters": {
            "g": {
                "A": {"active": True, "min": 5, "max": 1},                # inverted bounds
                "B": {"active": True, "min": 0, "max": 10, "start": 50},  # start out of bounds
                "C": {"active": True, "min": 0, "max": 1, "prior": {"dist": "cauchy"}},
            }
        },
    }
    with pytest.raises(ValueError) as ei:
        cfgmod.validate_config(bad)
    msg = str(ei.value)
    for needle in ("method.preset", "validation.scheme", "objective.weighting",
                   "execution.backend", "assimilation.mode", "gating.cultivar",
                   "num_cores", "min (5) must be < max (1)", "start (50) is outside",
                   "prior.dist 'cauchy'"):
        assert needle in msg, needle


def test_validate_config_requires_an_active_parameter():
    cfg = {"parameters": {"g": {"A": {"active": False, "min": 0, "max": 1}}}}
    with pytest.raises(ValueError, match="No active parameters"):
        cfgmod.validate_config(cfg)


@pytest.mark.parametrize(
    ("group", "level", "gate"),
    [
        ("genetic_cultivar", "cultivar", "blocked"),
        ("genetic_ecotype", "ecotype", "blocked"),
        ("genetic_species", "species", "blocked"),
        ("genetic_species", "species", "gated"),
    ],
)
def test_validate_config_rejects_active_parameters_behind_closed_gates(
    group, level, gate
):
    cfg = {
        "gating": {level: gate},
        "parameters": {
            group: {
                "P": {
                    "active": True,
                    "min": 0.0,
                    "max": 1.0,
                    "start": 0.5,
                }
            }
        },
    }
    with pytest.raises(ValueError, match=rf"gating\.{level}"):
        cfgmod.validate_config(cfg)


def test_load_config_validate_false_skips_check(tmp_path):
    # A config with no active parameters loads fine when validation is disabled.
    p = tmp_path / "c.yaml"
    p.write_text("parameters: {}\n")
    cfg = cfgmod.load_config(p, validate=False)
    assert cfg["parameters"] == {}


def test_load_hemp_config():
    cfg = cfgmod.load_config(HEMP_CFG)
    assert cfg["calibrator"]["name"] == "hemp_yunnan_2021_2022"
    assert cfg["method"]["preset"] == "A"
    crop = cfgmod.crop_for(cfg, "HM")
    assert crop["model"] == "CRGRO048"
    assert crop["cultivar_anchor"] == "IB0008"


def test_active_parameters():
    cfg = cfgmod.load_config(HEMP_CFG)
    active = cfgmod.active_parameters(cfg)
    names = {p["name"] for p in active}
    assert names == {
        "CSDL", "PPSEN", "EM-FL", "FL-SH", "FL-SD", "SD-PM", "FL-LF",
        "LFMAX", "SLAVR", "SIZLF", "XFRT", "WTPSD", "SFDUR", "SDPDV",
        "PODUR", "THRSH", "THVAR", "PL-EM", "EM-V1", "V1-JU", "JU-R0",
        "PM09", "LNGSH", "FL-VS", "TRIFL", "RWDTH", "RHGHT", "R1PPO",
        "OPTBI", "SLOBI",
    }
    for p in active:
        assert p["min"] < p["max"]
        assert p["min"] <= p["start"] <= p["max"], p["name"]


def test_cultivar_scoped_parameters_expand_to_configured_cultivars():
    cfg = deepcopy(cfgmod.load_config(HEMP_CFG))
    cfg["crops"][0]["calibration_cultivars"] = ["IB0008"]
    for group, params in cfg["parameters"].items():
        for spec in params.values():
            if isinstance(spec, dict):
                spec["active"] = False
    cfg["parameters"]["genetic_cultivar"]["EM-FL"].update(
        {"active": True, "scope": "cultivar"}
    )
    space = ParameterSpace.from_config(cfg)
    assert space.names == ["EM-FL__IB0008"]
    assert space.specs[0]["base_name"] == "EM-FL"
    assert space.specs[0]["cultivar"] == "IB0008"


def test_group_default_scope_expands_cultivar_parameters():
    cfg = deepcopy(cfgmod.load_config(HEMP_CFG))
    cfg["crops"][0]["calibration_cultivars"] = ["IB0008"]
    cfg["parameter_defaults"] = {"scope_by_group": {"genetic_cultivar": "cultivar"}}
    for group, params in cfg["parameters"].items():
        for spec in params.values():
            if isinstance(spec, dict):
                spec["active"] = False
                spec.pop("scope", None)
    cfg["parameters"]["genetic_cultivar"]["EM-FL"]["active"] = True
    space = ParameterSpace.from_config(cfg)
    assert space.names == ["EM-FL__IB0008"]


def test_cultivar_scoped_parameters_use_per_cultivar_starts():
    cfg = deepcopy(cfgmod.load_config(HEMP_CFG))
    cfg["crops"][0]["calibration_cultivars"] = ["IB0002", "IB0008"]
    for group, params in cfg["parameters"].items():
        for spec in params.values():
            if isinstance(spec, dict):
                spec["active"] = False
    cfg["parameters"]["genetic_cultivar"]["EM-FL"].update({
        "active": True,
        "scope": "cultivar",
        "max": 110.0,
        "start_by_cultivar": {"IB0002": 60.0, "IB0008": 85.0},
    })

    space = ParameterSpace.from_config(cfg)

    assert space.names == ["EM-FL__IB0002", "EM-FL__IB0008"]
    assert list(space.start) == [60.0, 85.0]


def test_active_cultivar_subset_can_keep_other_cultivars_fixed():
    cfg = deepcopy(cfgmod.load_config(HEMP_CFG))
    cfg["crops"][0]["calibration_cultivars"] = ["IB0002", "IB0008"]
    for group, params in cfg["parameters"].items():
        for spec in params.values():
            if isinstance(spec, dict):
                spec["active"] = False
                spec["fixed"] = False
    cfg["parameters"]["genetic_cultivar"]["EM-FL"].update({
        "active": True,
        "fixed": True,
        "scope": "cultivar",
        "cultivars": ["IB0008"],
        "fixed_cultivars": ["IB0002"],
        "start_by_cultivar": {"IB0002": 60.0, "IB0008": 85.0},
    })

    active = ParameterSpace.from_config(cfg)
    fixed = expand_parameter_specs(cfg, cfgmod.fixed_parameters(cfg))

    assert active.names == ["EM-FL__IB0008"]
    assert [spec["name"] for spec in fixed] == ["EM-FL__IB0002"]
    assert fixed[0]["start"] == 60.0
    assert fixed[0]["fixed"] is True


def test_validate_config_rejects_scoped_bounds_outside_declared_bounds():
    cfg = {
        "parameters": {
            "genetic_cultivar": {
                "P": {
                    "active": True,
                    "min": 0.0,
                    "max": 1.0,
                    "start": 0.5,
                    "scope": "cultivar",
                    "min_by_cultivar": {"C1": -0.1},
                    "max_by_cultivar": {"C1": 1.2},
                    "start_by_cultivar": {"C1": 0.6},
                }
            }
        },
        "crops": [{"calibration_cultivars": ["C1"]}],
    }

    with pytest.raises(ValueError, match="must stay within declared bounds"):
        cfgmod.validate_config(cfg)


def test_parameter_space_reports_values_at_declared_writer_step():
    cfg = {
        "parameters": {
            "genetic_species": {
                "TB": {
                    "active": True,
                    "min": 1.47,
                    "max": 1.56,
                    "start": 1.514,
                    "step": 0.1,
                }
            }
        },
        "gating": {"species": "free"},
    }

    space = ParameterSpace.from_config(cfg)

    assert space.start.tolist() == pytest.approx([1.5])
    assert space.to_theta([1.556])["TB"] == pytest.approx(1.5)


def test_parameter_space_rejects_range_without_writer_grid_value():
    cfg = {
        "parameters": {
            "genetic_species": {
                "TB": {
                    "active": True,
                    "min": 1.51,
                    "max": 1.59,
                    "start": 1.55,
                    "step": 0.1,
                }
            }
        },
        "gating": {"species": "free"},
    }

    with pytest.raises(ValueError, match="contains no value"):
        ParameterSpace.from_config(cfg)


def test_grid_sampler_can_skip_start_row_for_exact_factorial_size():
    cfg = deepcopy(cfgmod.load_config(HEMP_CFG))
    cfg["crops"][0]["calibration_cultivars"] = ["IB0008"]
    for group, params in cfg["parameters"].items():
        for spec in params.values():
            if isinstance(spec, dict):
                spec["active"] = False
    for name in ("CSDL", "PPSEN"):
        cfg["parameters"]["genetic_cultivar"][name].update(
            {"active": True, "scope": "cultivar"}
        )
    space = ParameterSpace.from_config(cfg)
    design = sample(space, n=25, engine="grid", include_start=False)
    assert len(design) == 25
    assert design["CSDL__IB0008"].nunique() == 5
    assert design["PPSEN__IB0008"].nunique() == 5


def test_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv("DSSATCAL_NUM_CORES", "5")
    monkeypatch.setenv("DSSATCAL_HEMP_DIR", str(tmp_path))
    cfg = cfgmod.load_config(HEMP_CFG)
    assert cfg["calibrator"]["num_cores"] == 5
    assert cfg["source"]["hemp_dir"] == str(tmp_path)


def test_zero_num_cores_uses_all_logical_cores(monkeypatch):
    monkeypatch.setattr("os.cpu_count", lambda: 16)
    assert resolve_cores(0) == 16


def test_resolve_exe():
    cfg = cfgmod.load_config(HEMP_CFG)
    exe = cfgmod.resolve_exe(cfg)
    # The portable example has no machine-specific executable. Resolve a
    # discovered sibling install when available, otherwise return the canonical
    # DSSAT executable path that validation can report to the user.
    assert exe.name.lower() in {"dscsm048", "dscsm048.exe"}


def test_resolve_exe_retains_custom_name_when_native_discovery_fails(monkeypatch):
    custom = Path(
        "C:/Users/example/DSSAT48Hemp/"
        "dscsm048_compiled_4.8.2.with_HM_code.exe"
    )
    cfg = {
        "calibrator": {
            "dssat_dir": "C:/Users/example/DSSAT48Hemp",
            "dssat_exe": str(custom),
        }
    }
    monkeypatch.setattr(cfgmod, "_workspace_dssat_root", lambda: None)

    assert cfgmod.resolve_exe(cfg) == custom


def test_shared_stack_defaults_and_template_env(monkeypatch, tmp_path):
    assert cfgmod.DEFAULTS["execution"]["backend"] == "native"
    assert cfgmod.DEFAULTS["soil"]["provider"] == "file"

    monkeypatch.setenv("DSSAT_TEMPLATE_DIR", str(tmp_path))
    assert cfgmod.resolve_template_dir({}) == tmp_path

    cfg = cfgmod.load_config(HEMP_CFG)
    paths = cfgmod.resolve_dssat_paths(cfg)
    assert paths["genotype"].name == "Genotype"
    assert paths["weather"].name == "Weather"
    assert paths["soil"].name == "Soil"
