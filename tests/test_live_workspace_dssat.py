"""Portable real-DSSAT integration check against the sibling workspace install."""

from pathlib import Path

import pytest

from dssatcalibrator.config import active_parameters, crop_for, load_config, resolve_exe
from dssatcalibrator.spawn import spawn_and_run


@pytest.mark.slow
def test_current_workspace_hemp_spawn_selected_treatment(tmp_path):
    repo = Path(__file__).resolve().parents[1]
    dssat_root = repo.parent / "DSSAT48"
    hemp_dir = dssat_root / "Hemp"
    experiment = "UFCI2101"
    required = [
        dssat_root / "dscsm048",
        hemp_dir / f"{experiment}.HMX",
        dssat_root / "Genotype" / "HMGRO048.CUL",
    ]
    if not all(path.exists() for path in required):
        pytest.skip("sibling DSSAT48 hemp integration fixture is not installed")

    cfg = load_config(repo / "config_hemp.yaml", validate=False)
    cfg["calibrator"]["dssat_dir"] = str(dssat_root)
    cfg["calibrator"]["dssat_exe"] = str(dssat_root / "dscsm048")
    cfg["source"]["hemp_dir"] = str(hemp_dir)
    cfg["experiments"] = [experiment]
    crop = crop_for(cfg, "HM")
    crop.update({
        "cultivar_anchor": "IB0001",
        "cultivar_anchors": ["IB0001", "IB0002"],
        "ecotype": "HM0001",
        "cultivar_ecotypes": {"IB0001": "HM0001", "IB0002": "HM0002"},
    })
    specs = active_parameters(cfg)
    theta = {spec["name"]: spec["start"] for spec in specs}
    result = spawn_and_run(
        theta, exp_id=experiment, cfg=cfg, crop=crop, param_specs=specs,
        run_root=tmp_path, treatments=[2], exe=resolve_exe(cfg),
    )

    assert result.status == "success", result.message
    assert not result.plantgro.empty and (result.plantgro["CWAD"] > 0).any()
    assert set(result.plantgro["treatment"].dropna().astype(int)) == {2}
    assert set(result.evaluate["treatment"].dropna().astype(int)) == {2}
