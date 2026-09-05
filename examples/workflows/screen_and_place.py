#!/usr/bin/env python3
"""Screen metadata, resolve candidate CIFs, place CO2, and retain per-candidate outcomes.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, stage_demo_data, write_report


def main():
    parser = example_parser("screen_and_place", __doc__, seed=True)
    parser.add_argument("--data-path", type=Path)
    parser.add_argument("--structures-dir", type=Path)
    parser.add_argument("--metal", default="Zn")
    parser.add_argument("--limit", type=int, default=3)
    args = parser.parse_args()
    if bool(args.data_path) != bool(args.structures_dir):
        parser.error("Real mode requires --data-path and --structures-dir together.")
    out = output_dir(args)
    import shutil

    from mofforge import Crystal, infer_bonds, place_adsorbate
    from mofforge.coremof import CoreMOFDatabase
    from mofforge.coremof.structures import resolve_structure_path

    if args.data_path:
        source = out / "coremof_input.csv"
        shutil.copyfile(args.data_path, source)
        structures = args.structures_dir
    else:
        _, source, structures = stage_demo_data(out)
    db = CoreMOFDatabase(source)
    try:
        candidates = db.screen(
            metal=args.metal,
            pld_min=3.8,
            water_stability_min=0.7,
            extension="All Solvent Removed",
            limit=args.limit,
        )
    finally:
        db.close()
    results = []
    for rec in candidates:
        source_cif = resolve_structure_path(rec.coreid, structures_dir=structures)
        if source_cif is None:
            results.append({"coreid": rec.coreid, "status": "missing_structure"})
            continue
        try:
            host = infer_bonds(Crystal.from_cif(source_cif))
            placed = place_adsorbate(
                host, "CO2", n_adsorbates=2, grid_spacing=1.0, random_seed=args.seed
            )
            destination = out / f"{rec.coreid}_CO2.cif"
            placed.crystal.write_cif(destination)
            results.append(
                {
                    "coreid": rec.coreid,
                    "status": "placed",
                    "output": destination,
                    "requested": 2,
                    "placed": placed.n_adsorbates,
                    "validation": placed.validation.to_dict(),
                }
            )
        except (OSError, ValueError) as exc:
            results.append({"coreid": rec.coreid, "status": "failed", "error": str(exc)})
    write_report(
        out,
        synthetic_data=args.data_path is None,
        seed=args.seed,
        candidates=len(candidates),
        results=results,
    )
    if not any(row["status"] == "placed" for row in results):
        raise SystemExit("No candidate produced a loaded structure; inspect summary.json.")


if __name__ == "__main__":
    main()
