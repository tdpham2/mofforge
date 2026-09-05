#!/usr/bin/env python3
"""Replace selected phenylene hydrogens with acetylamido groups.

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
    parser = example_parser("linker_functionalization", __doc__, seed=True)
    parser.add_argument("--crystal", type=Path, default=CRYSTALS / "IRMOF-1.cif")
    parser.add_argument("--nb-loc", type=int, default=6)
    parser.add_argument("--output", help="Override the primary CIF filename.")
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

    parent = infer_bonds(Crystal.from_cif(args.crystal))
    # H! marks the atom to remove; untagged atoms provide the alignment scaffold.
    query = fragment("2-!-p-phenylene.xyz", fragment_path=args.fragment_path)
    replacement = fragment("2-acetylamido-p-phenylene.xyz", fragment_path=args.fragment_path)
    match = find_pattern(query, parent)
    if not match.nb_locations():
        raise SystemExit("Query not found: check explicit linker hydrogens.")
    child = replace_pattern(match, replacement, nb_loc=args.nb_loc, random_seed=args.seed)
    destination = output_file(args, "acetylamido_IRMOF-1.cif")
    child.write_cif(destination)
    write_report(
        out,
        seed=args.seed,
        parent_atoms=parent.n_atoms,
        child_atoms=child.n_atoms,
        locations=child.provenance.parameters["locations"],
        output=destination,
        validation=validate_structure(infer_bonds(child)).to_dict(),
    )


if __name__ == "__main__":
    main()
