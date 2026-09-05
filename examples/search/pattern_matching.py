#!/usr/bin/env python3
"""Find phenylene locations and interpret atom mappings.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, MOIETIES, example_parser, output_dir, write_report


def main():
    parser = example_parser("pattern_matching", __doc__)
    parser.add_argument("--crystal", type=Path, default=CRYSTALS / "IRMOF-1.cif")
    parser.add_argument("--query", default="p-phenylene.xyz")
    parser.add_argument("--fragment-path", type=Path, default=MOIETIES)
    args = parser.parse_args()
    out = output_dir(args)
    from mofforge import Crystal, find_pattern, fragment, infer_bonds

    parent = infer_bonds(Crystal.from_cif(args.crystal), periodic=True)
    query = fragment(args.query, fragment_path=args.fragment_path)
    # Locations are atom sets; orientations are mappings within each set.
    match = find_pattern(query, parent)
    if not match.nb_locations():
        raise SystemExit("No matches: check the query, hydrogens, and inferred bonds.")
    match.matched_substructures().write_xyz(out / "matched_atoms.xyz")
    write_report(
        out,
        atoms=parent.n_atoms,
        locations=match.nb_locations(),
        isomorphisms=match.nb_isomorphisms(),
        orientations=match.nb_ori_at_loc(),
        first_mapping=match.isomorphisms[0][0],
    )


if __name__ == "__main__":
    main()
