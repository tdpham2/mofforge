#!/usr/bin/env python3
"""Compare all five replacement modes on the same parent.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, MOIETIES, example_parser, output_dir, write_report


def main():
    args = example_parser("selective_modification", __doc__, seed=True).parse_args()
    out = output_dir(args)
    from mofforge import (
        Crystal,
        find_pattern,
        fragment,
        infer_bonds,
        replace_pattern,
        validate_structure,
    )

    parent = infer_bonds(Crystal.from_cif(CRYSTALS / "IRMOF-1.cif"))
    query = fragment("2-!-p-phenylene.xyz", fragment_path=MOIETIES)
    replacement = fragment("2-nitro-p-phenylene.xyz", fragment_path=MOIETIES)
    match = find_pattern(query, parent)
    # Both loc and explicit ori indices are zero-based. Omit ori to optimize.
    modes = {
        "all_optimal": {},
        "random_locations": {"nb_loc": 8},
        "specific_locations": {"loc": [0, 5, 10, 15]},
        "specific_orientations": {"loc": [0, 1, 2, 3], "ori": [0, 1, 2, 3]},
        "random_orientations": {"random": True},
    }
    rows = {}
    for name, options in modes.items():
        child = replace_pattern(match, replacement, name=name, random_seed=args.seed, **options)
        child.write_cif(out / f"{name}.cif")
        rows[name] = {
            "atoms": child.n_atoms,
            "locations": child.provenance.parameters["locations"],
            "orientations": child.provenance.parameters["orientations"],
            "validation": validate_structure(infer_bonds(child)).to_dict(),
        }
    write_report(out, seed=args.seed, parent_atoms=parent.n_atoms, modes=rows)


if __name__ == "__main__":
    main()
