"""Run a real-DSSAT one-at-a-time parameter impact atlas.

Example:

    python run_impact_atlas.py config_hemp.yaml --experiments UFCI2101 UFCI2201 \
        --groups genetic_cultivar genetic_ecotype genetic_species --discover-genotype \
        --allow-species \
        --max-parameters 12 --cores 4 --no-long
"""

from __future__ import annotations

import argparse
from pathlib import Path

from dssatcalibrator.config import load_config
from dssatcalibrator.impact import DEFAULT_GROUPS, run_impact_atlas


def main() -> None:
    ap = argparse.ArgumentParser(description="Run a real-DSSAT parameter impact atlas.")
    ap.add_argument("config", help="calibration YAML")
    ap.add_argument("--outdir", default=None, help="output directory")
    ap.add_argument("--experiments", nargs="*", default=None, help="experiment IDs to include")
    ap.add_argument("--groups", nargs="*", default=None,
                    help=f"parameter groups to sweep (default: {' '.join(DEFAULT_GROUPS)})")
    ap.add_argument("--levels", nargs="*", default=["low", "high"],
                    help="levels to run for each parameter: low high start")
    ap.add_argument("--active-only", action="store_true",
                    help="only sweep parameters with active: true in the config")
    ap.add_argument("--discover-cultivar", action="store_true",
                    help="discover cultivar coefficients from the crop .CUL file")
    ap.add_argument("--discover-ecotype", action="store_true",
                    help="discover ecotype coefficients from the crop .ECO file")
    ap.add_argument("--discover-species", action="store_true",
                    help="discover numeric species tokens from the crop .SPE file")
    ap.add_argument("--discover-genotype", action="store_true",
                    help="discover cultivar, ecotype, and species genotype parameters")
    ap.add_argument("--allow-species", action="store_true",
                    help="explicitly allow .SPE species edits for species sweeps")
    ap.add_argument("--max-parameters", type=int, default=None,
                    help="limit number of parameters for smoke testing")
    ap.add_argument("--max-per-group", type=int, default=None,
                    help="limit number of parameters per group for smoke testing")
    ap.add_argument("--outputs", nargs="*", default=None,
                    help="specific DSSAT *.OUT files to collect")
    ap.add_argument("--no-long", action="store_true",
                    help="do not write outputs_long.csv; effects are still computed")
    ap.add_argument("--compress-long", action="store_true",
                    help="write outputs_long.csv.gz instead of outputs_long.csv")
    ap.add_argument("--effect-tolerance", type=float, default=1e-9,
                    help="absolute delta threshold used to count changed outputs")
    ap.add_argument("--cores", type=int, default=None,
                    help="override calibrator.num_cores")
    ap.add_argument("--dssat-exe", default=None, help="override calibrator.dssat_exe")
    ap.add_argument("--dssat-dir", default=None, help="override calibrator.dssat_dir")
    ap.add_argument("--hemp-dir", default=None, help="override source.hemp_dir")
    ap.add_argument("--keep-existing", action="store_true",
                    help="do not delete the output directory before running")
    ap.add_argument("--no-progress", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config, validate=False)
    if args.dssat_exe:
        cfg.setdefault("calibrator", {})["dssat_exe"] = args.dssat_exe
    if args.dssat_dir:
        cfg.setdefault("calibrator", {})["dssat_dir"] = args.dssat_dir
    if args.hemp_dir:
        cfg.setdefault("source", {})["hemp_dir"] = args.hemp_dir
    name = cfg.get("calibrator", {}).get("name", "run")
    outdir = Path(args.outdir or Path(cfg.get("calibrator", {}).get("results_dir", "results")) / f"{name}_impact_atlas")

    result = run_impact_atlas(
        cfg,
        output_dir=outdir,
        groups=args.groups,
        experiments=args.experiments,
        levels=args.levels,
        active_only=args.active_only,
        discover_cultivar=args.discover_cultivar,
        discover_ecotype=args.discover_ecotype,
        discover_species=args.discover_species,
        discover_genotype=args.discover_genotype,
        allow_species=args.allow_species,
        max_parameters=args.max_parameters,
        max_per_group=args.max_per_group,
        output_files=args.outputs,
        num_cores=args.cores,
        write_long=not args.no_long,
        compress_long=args.compress_long,
        effect_tolerance=args.effect_tolerance,
        clean=not args.keep_existing,
        progress=not args.no_progress,
    )

    ok = result.run_manifest["status"].isin(["success", "cached"]).sum()
    total = len(result.run_manifest)
    print(f"\nImpact atlas complete: {ok}/{total} spawns succeeded")
    print(f"Output directory: {result.output_dir.resolve()}")
    print(f"Run manifest: {result.output_dir / 'run_manifest.csv'}")
    print(f"Parameter catalog: {result.output_dir / 'parameter_catalog.csv'}")
    print(f"Effects table: {result.output_dir / 'parameter_output_effects.csv'}")
    print(f"Impact summary: {result.output_dir / 'impact_summary.md'}")
    print(f"Parameter summary: {result.output_dir / 'parameter_impact_summary.csv'}")
    print(f"Capability map: {result.output_dir / 'capability_map.md'}")
    if not args.no_long:
        suffix = "outputs_long.csv.gz" if args.compress_long else "outputs_long.csv"
        print(f"Long output: {result.output_dir / suffix}")


if __name__ == "__main__":
    main()
