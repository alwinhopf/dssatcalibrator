"""CLI: scaffold a new crop/cultivar from an analog DSSAT module.

Clones the analog genotype files under a new stem/code and prints (and optionally
writes) a starter ``parameters:`` config block you then review against literature.

    python scaffold_crop.py --dssat-dir C:/DSSAT48 \
        --analog-stem SBGRO048 --new-stem QUGRO048 --new-code QU \
        --source-anchor IB0001 --out-dir templates/quinoa

See CONCEPT.md §18 / WALKTHROUGH.md §17 for the new-crop workflow.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from dssatcalibrator.scaffold import scaffold_crop


def main():
    ap = argparse.ArgumentParser(description="Scaffold a new crop from an analog DSSAT module.")
    ap.add_argument("--dssat-dir", required=True, help="DSSAT install dir (contains Genotype/)")
    ap.add_argument("--analog-stem", required=True, help="analog genotype stem, e.g. SBGRO048")
    ap.add_argument("--new-stem", required=True, help="new genotype stem, e.g. QUGRO048")
    ap.add_argument("--new-code", required=True, help="new 2-letter crop code, e.g. QU")
    ap.add_argument("--source-anchor", required=True, help="cultivar row to copy starts from, e.g. IB0001")
    ap.add_argument("--new-anchor", default=None, help="new cultivar anchor code (default: source)")
    ap.add_argument("--out-dir", default=None,
                    help="where to write the cloned .CUL/.ECO/.SPE "
                         "(default: DSSAT_TEMPLATE_DIR or shared dssat_templates)")
    ap.add_argument("--spread", type=float, default=0.3,
                    help="+/- fraction for bounds when no MINIMA/MAXIMA rows exist (default 0.3)")
    ap.add_argument("--no-spe", action="store_true", help="do not copy the .SPE")
    args = ap.parse_args()

    res = scaffold_crop(
        dssat_dir=args.dssat_dir, analog_stem=args.analog_stem, new_stem=args.new_stem,
        new_code=args.new_code, source_anchor=args.source_anchor, new_anchor=args.new_anchor,
        out_dir=args.out_dir, spread=args.spread, copy_spe=not args.no_spe,
    )

    print("Cloned genotype files:")
    for ext, p in res["files"].items():
        print(f"  .{ext}: {p}")
    block_path = Path(res["out_dir"]) / "parameters_block.yaml"
    block_path.write_text(res["parameters_yaml"], encoding="utf-8")
    print(f"\nStarter parameter block -> {block_path}")
    print("Paste it into your config and review every bound against literature.\n")
    print(res["parameters_yaml"])


if __name__ == "__main__":
    main()
