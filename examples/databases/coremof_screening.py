#!/usr/bin/env python3
"""Search CoRE metadata and screen by pore dimensions and stability.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, stage_demo_data, write_report


def main():
    parser = example_parser("coremof_screening", __doc__)
    parser.add_argument("--data-path", type=Path)
    parser.add_argument("--metal", default="Zn")
    parser.add_argument("--pld-min", type=float, default=3.8)
    parser.add_argument("--water-stability-min", type=float, default=0.7)
    args = parser.parse_args()
    out = output_dir(args)
    import shutil
    from dataclasses import asdict

    from mofforge.coremof import CoreMOFDatabase

    if args.data_path:
        source = out / "coremof_input.csv"
        shutil.copyfile(args.data_path, source)
    else:
        _, source, _ = stage_demo_data(out)
    db = CoreMOFDatabase(source)
    try:
        search = db.search(args.metal, field="metal", limit=5)
        candidates = db.screen(
            metal=args.metal,
            pld_min=args.pld_min,
            water_stability_min=args.water_stability_min,
            extension="All Solvent Removed",
            limit=5,
        )
        empty = db.screen(pld_min=10000, limit=5)
        write_report(
            out,
            synthetic_data=args.data_path is None,
            records=[asdict(r) for r in search.records],
            candidates=[asdict(r) for r in candidates],
            deliberately_empty_count=len(empty),
            cached_metadata=source,
        )
    finally:
        db.close()


if __name__ == "__main__":
    main()
