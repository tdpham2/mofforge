#!/usr/bin/env python3
"""Bridge CSD names to CoRE variants and resolve separately stored CIF files.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, stage_demo_data, write_report


def main():
    parser = example_parser("bridge_and_resolve", __doc__)
    parser.add_argument("--csd-data", type=Path)
    parser.add_argument("--coremof-data", type=Path)
    parser.add_argument("--structures-dir", type=Path)
    parser.add_argument("--name", default="Fictional")
    args = parser.parse_args()
    if any([args.csd_data, args.coremof_data, args.structures_dir]) and not all(
        [args.csd_data, args.coremof_data, args.structures_dir]
    ):
        parser.error(
            "Real mode requires --csd-data, --coremof-data, and --structures-dir together."
        )
    out = output_dir(args)
    import shutil

    from mofforge.coremof import CoreMOFDatabase, search_csd_name
    from mofforge.coremof.structures import resolve_structure_path
    from mofforge.csd import CSDDatabase

    if args.csd_data:
        csd_path, core_path = out / "csd_input.tab", out / "coremof_input.csv"
        shutil.copyfile(args.csd_data, csd_path)
        shutil.copyfile(args.coremof_data, core_path)
        structures = args.structures_dir
    else:
        csd_path, core_path, structures = stage_demo_data(out)
    csd, core = CSDDatabase(csd_path), CoreMOFDatabase(core_path)
    try:
        matches = search_csd_name(args.name, csd_db=csd, coremof_db=core, limit=5)
        rows = []
        for match in matches:
            variants = []
            for rec in match.coremof_records:
                cif = resolve_structure_path(rec.coreid, structures_dir=structures)
                variants.append({"coreid": rec.coreid, "extension": rec.extension, "cif": cif})
            rows.append({"refcode": match.csd_record.refcode, "variants": variants})
        write_report(out, synthetic_data=args.csd_data is None, matches=rows)
    finally:
        csd.close()
        core.close()


if __name__ == "__main__":
    main()
