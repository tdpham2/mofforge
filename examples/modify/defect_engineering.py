#!/usr/bin/env python3
"""Remove UiO-66 linker cores while installing formate caps.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import (
    CRYSTALS,
    MOIETIES,
    example_parser,
    output_dir,
    output_file,
    write_report,
)


def main():
    parser = example_parser("defect_engineering", __doc__, seed=True)
    parser.add_argument("--loc", nargs="+", type=int, default=[2, 8])
    parser.add_argument("--output")
    parser.add_argument("--fragment-path", type=Path, default=MOIETIES)
    args = parser.parse_args()
    out = output_dir(args)
    from mofforge import (
        Crystal,
        find_pattern,
        fragment,
        infer_bonds,
        replace_pattern,
        validate_structure,
    )

    parent = infer_bonds(Crystal.from_cif(CRYSTALS / "UiO-66.cif"))
    query = fragment("BDC.xyz", fragment_path=args.fragment_path)
    caps = fragment("formate_caps.xyz", fragment_path=args.fragment_path)
    match = find_pattern(query, parent)
    if not match.nb_locations():
        raise SystemExit("Expected BDC linkers were not found.")
    child = replace_pattern(match, caps, loc=args.loc, random_seed=args.seed)
    destination = output_file(args, "defected_UiO-66.cif")
    child.write_cif(destination)
    write_report(
        out,
        locations=args.loc,
        matches=match.nb_locations(),
        parent_atoms=parent.n_atoms,
        child_atoms=child.n_atoms,
        output=destination,
        validation=validate_structure(infer_bonds(child)).to_dict(),
    )


if __name__ == "__main__":
    main()
