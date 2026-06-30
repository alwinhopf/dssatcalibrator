from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from dssatcalibrator.dssat_io import (
    collect_run_outputs,
    dssat_output_long,
    parse_dssat_output_tables,
)
from dssatcalibrator import impact
from dssatcalibrator.impact import (
    capability_map,
    discover_cultivar_parameters,
    discover_ecotype_parameters,
    discover_species_parameters,
)


def test_generic_dssat_output_parser_to_long(tmp_path):
    out = tmp_path / "SoilWat.OUT"
    out.write_text(
        "*SOIL WATER\n"
        "@YEAR DOY   DAS   SWTD   ESAD\n"
        " 2021 120     1   123.4   1.2\n"
        " 2021 121     2   -99     1.4\n",
        encoding="utf-8",
    )

    wide = parse_dssat_output_tables(out)
    assert wide["SWTD"].iloc[0] == 123.4
    assert pd.isna(wide["SWTD"].iloc[1])
    assert wide["date"].dt.strftime("%Y-%m-%d").tolist() == ["2021-04-30", "2021-05-01"]

    long = dssat_output_long(wide)
    assert {"SWTD", "ESAD"} <= set(long["variable"])
    assert long.query("variable == 'SWTD'")["value"].tolist() == [123.4]


def test_generic_dssat_output_parser_preserves_text_leading_ids(tmp_path):
    out = tmp_path / "PlantNBal.OUT"
    out.write_text(
        "*N BALANCE\n"
        "@EXCODE TRNO MODEL CWAD NUPC\n"
        " UFCI2101   1 CRGRO 123.4  5.6\n",
        encoding="utf-8",
    )

    wide = parse_dssat_output_tables(out)
    assert wide["EXCODE"].tolist() == ["UFCI2101"]

    long = dssat_output_long(wide)
    assert {"excode", "model"} <= set(long.columns)
    assert long.query("variable == 'CWAD'")["value"].tolist() == [123.4]
    assert long.query("variable == 'CWAD'")["excode"].tolist() == ["UFCI2101"]

    collected = collect_run_outputs(tmp_path)
    manifest = dict(zip(collected["manifest"]["source_file"], collected["manifest"]["exists"]))
    assert bool(manifest["PlantNBal.OUT"]) is True


def test_capability_map_routes_reusable_io_to_support_packages():
    cap = capability_map([
        {"group": "genetic_cultivar", "name": "CSDL"},
        {"group": "soil", "name": "SDUL"},
        {"group": "weather", "name": "SRAD"},
    ])
    owners = dict(zip(cap["parameter"], cap["recommended_owner"]))
    assert owners["CSDL"] == "dssatengine"
    assert owners["SDUL"] == "dssatengine + dssatutils"
    assert owners["SRAD"] == "dssatengine + dssatutils"


def test_discover_species_parameters_from_spe(tmp_path, monkeypatch):
    genotype = tmp_path / "Genotype"
    genotype.mkdir()
    (genotype / "HMGRO048.SPE").write_text(
        "  40.00 61.00  0.96  0.10                   PARMAX,PHTMAX,KCAN,KC_SLOPE\n",
        encoding="utf-8",
    )
    cfg = {
        "calibrator": {"dssat_dir": str(tmp_path)},
        "crops": [{"code": "HM", "genotype_stem": "HMGRO048"}],
    }
    specs = discover_species_parameters(cfg)
    assert [s["name"] for s in specs] == ["PARMAX", "PHTMAX", "KCAN", "KC_SLOPE"]
    assert specs[2]["spe_index"] == 2
    assert specs[2]["spe_key"] == "PARMAX,PHTMAX,KCAN,KC_SLOPE"


def test_discover_cultivar_and_ecotype_parameters_from_genotype_files(tmp_path):
    genotype = tmp_path / "Genotype"
    genotype.mkdir()
    fixture_dir = Path("tests/fixtures")
    (genotype / "HMGRO048.CUL").write_text((fixture_dir / "sample.CUL").read_text(), encoding="utf-8")
    (genotype / "HMGRO048.ECO").write_text((fixture_dir / "sample.ECO").read_text(), encoding="utf-8")
    cfg = {
        "calibrator": {"dssat_dir": str(tmp_path)},
        "crops": [{
            "code": "HM",
            "genotype_stem": "HMGRO048",
            "cultivar_anchor": "IB0008",
            "ecotype": "SB0301",
        }],
    }

    cultivar = discover_cultivar_parameters(cfg)
    ecotype = discover_ecotype_parameters(cfg)

    by_name = {s["name"]: s for s in cultivar}
    assert by_name["CSDL"]["start"] == 12.33
    assert by_name["CSDL"]["min"] == 11.0
    assert by_name["CSDL"]["max"] == 14.0
    assert by_name["LFMAX"]["source"] == "discovered_cul"
    assert {s["name"] for s in ecotype} >= {"PL-EM", "JU-R0", "PM06"}
    assert all(s["source"] == "discovered_eco" for s in ecotype)


def test_run_impact_atlas_can_skip_long_output_and_write_catalog(tmp_path, monkeypatch):
    class FakeObs:
        table = pd.DataFrame()

        def experiments(self):
            return ["EXP1"]

    def fake_run_many(jobs, n_workers, on_done=None):
        out = []
        for i, job in enumerate(jobs):
            res = SimpleNamespace(
                status="success",
                run_dir=tmp_path / f"run_{i}",
                theta=job["theta"],
                message="",
            )
            if on_done:
                on_done(res)
            out.append(res)
        return out

    def fake_collect(_run_dir, output_files=None):
        long = pd.DataFrame({
            "source_file": ["PlantGro.OUT"],
            "variable": ["CWAD"],
            "treatment": [1],
            "row_index": [0],
            "value": [100.0],
        })
        manifest = pd.DataFrame({"source_file": ["PlantGro.OUT"], "exists": [True], "size_bytes": [10]})
        return {
            "wide": pd.DataFrame(),
            "long": long,
            "plantgro": pd.DataFrame(),
            "evaluate": pd.DataFrame(),
            "summary": pd.DataFrame(),
            "manifest": manifest,
        }

    monkeypatch.setattr(impact, "_load_observations", lambda cfg, experiments, crop: FakeObs())
    monkeypatch.setattr(impact, "parse_treatments", lambda path: [1])
    monkeypatch.setattr(impact, "run_many", fake_run_many)
    monkeypatch.setattr(impact.dssat_io, "collect_run_outputs", fake_collect)
    monkeypatch.setattr(impact.obj, "score", lambda results, obs, cfg: SimpleNamespace(score=1.23, residuals=[1, 2]))

    cfg = {
        "calibrator": {"dssat_dir": str(tmp_path), "dssat_exe": str(tmp_path / "dssat.exe"), "num_cores": 1},
        "source": {"hemp_dir": str(tmp_path)},
        "crops": [{"code": "HM", "filex_ext": "HMX", "genotype_stem": "HMGRO048"}],
        "experiments": ["EXP1"],
        "parameters": {
            "weather": {
                "radiation_mult": {"active": False, "min": 0.9, "max": 1.1, "start": 1.0}
            }
        },
    }

    result = impact.run_impact_atlas(
        cfg,
        output_dir=tmp_path / "atlas",
        groups=["weather"],
        write_long=False,
        progress=False,
    )

    assert (tmp_path / "atlas" / "parameter_catalog.csv").exists()
    assert (tmp_path / "atlas" / "file_manifest.csv").exists()
    assert (tmp_path / "atlas" / "score_effects.csv").exists()
    assert (tmp_path / "atlas" / "output_impact_summary.csv").exists()
    assert (tmp_path / "atlas" / "parameter_impact_summary.csv").exists()
    assert (tmp_path / "atlas" / "impact_summary.md").exists()
    assert not (tmp_path / "atlas" / "outputs_long.csv").exists()
    assert result.output_long.empty
    assert not result.score_effects.empty
    assert not result.parameter_impact_summary.empty
    catalog = pd.read_csv(tmp_path / "atlas" / "parameter_catalog.csv")
    assert catalog["name"].tolist() == ["radiation_mult"]


def test_run_impact_atlas_errors_when_no_experiments_survive_filtering(tmp_path, monkeypatch):
    class FakeObs:
        table = pd.DataFrame()

        def experiments(self):
            return ["OTHER"]

    monkeypatch.setattr(impact, "_load_observations", lambda cfg, experiments, crop: FakeObs())
    cfg = {
        "calibrator": {"dssat_dir": str(tmp_path), "dssat_exe": str(tmp_path / "dssat.exe"), "num_cores": 1},
        "source": {"hemp_dir": str(tmp_path)},
        "crops": [{"code": "HM", "filex_ext": "HMX", "genotype_stem": "HMGRO048"}],
        "experiments": ["EXP1"],
        "parameters": {},
    }

    with pytest.raises(ValueError, match="no experiments"):
        impact.run_impact_atlas(cfg, output_dir=tmp_path / "atlas_empty", progress=False)


def test_run_impact_atlas_requires_explicit_species_permission(tmp_path, monkeypatch):
    class FakeObs:
        table = pd.DataFrame()

        def experiments(self):
            return ["EXP1"]

    seen_gates = []

    def fake_run_many(jobs, n_workers, on_done=None):
        out = []
        for i, job in enumerate(jobs):
            seen_gates.append(job["cfg"].get("gating", {}).get("species"))
            res = SimpleNamespace(
                status="success",
                run_dir=tmp_path / f"run_{i}",
                theta=job["theta"],
                message="",
            )
            if on_done:
                on_done(res)
            out.append(res)
        return out

    def fake_collect(_run_dir, output_files=None):
        return {
            "wide": pd.DataFrame(),
            "long": pd.DataFrame({
                "source_file": ["PlantGro.OUT"],
                "variable": ["CWAD"],
                "treatment": [1],
                "row_index": [0],
                "value": [100.0],
            }),
            "plantgro": pd.DataFrame(),
            "evaluate": pd.DataFrame(),
            "summary": pd.DataFrame(),
            "manifest": pd.DataFrame({"source_file": ["PlantGro.OUT"], "exists": [True], "size_bytes": [10]}),
        }

    monkeypatch.setattr(impact, "_load_observations", lambda cfg, experiments, crop: FakeObs())
    monkeypatch.setattr(impact, "parse_treatments", lambda path: [1])
    monkeypatch.setattr(impact, "run_many", fake_run_many)
    monkeypatch.setattr(impact.dssat_io, "collect_run_outputs", fake_collect)
    monkeypatch.setattr(impact.obj, "score", lambda results, obs, cfg: SimpleNamespace(score=1.0, residuals=[]))

    cfg = {
        "calibrator": {"dssat_dir": str(tmp_path), "dssat_exe": str(tmp_path / "dssat.exe"), "num_cores": 1},
        "source": {"hemp_dir": str(tmp_path)},
        "crops": [{"code": "HM", "filex_ext": "HMX", "genotype_stem": "HMGRO048"}],
        "experiments": ["EXP1"],
        "gating": {"species": "blocked"},
        "parameters": {
            "genetic_species": {
                "PARMAX": {"active": True, "min": 36.0, "max": 44.0, "start": 40.0, "spe_key": "PARMAX"}
            }
        },
    }

    with pytest.raises(ValueError, match="gating.species"):
        impact.run_impact_atlas(
            cfg,
            output_dir=tmp_path / "atlas_species_blocked",
            groups=["genetic_species"],
            write_long=False,
            progress=False,
        )

    impact.run_impact_atlas(
        cfg,
        output_dir=tmp_path / "atlas_species_allowed",
        groups=["genetic_species"],
        allow_species=True,
        write_long=False,
        progress=False,
    )

    assert seen_gates
    assert set(seen_gates) == {"free"}
