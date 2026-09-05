#!/usr/bin/env python3
"""Correct disordered rings and delete isolated acetylene guests.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, MOIETIES, example_parser, output_dir, write_report


def main():
    args = example_parser("cleanup_and_repair", __doc__, seed=True).parse_args()
    out = output_dir(args)
    from mofforge import (
        Crystal,
        find_pattern,
        fragment,
        infer_bonds,
        replace_pattern,
        validate_structure,
    )

    parent = infer_bonds(Crystal.from_cif(CRYSTALS / "SIFSIX-2-Cu-i.cif"))
    query = fragment("disordered_ligand!.xyz", fragment_path=MOIETIES)
    ring = fragment("4-pyridyl.xyz", fragment_path=MOIETIES)
    disorder = find_pattern(query, parent)
    if not disorder.nb_locations():
        raise SystemExit("Expected disordered rings were not found.")
    repaired = infer_bonds(replace_pattern(disorder, ring, random_seed=args.seed))
    repaired.write_cif(out / "ordered_with_guests.cif")
    # Match complete isolated molecules, protecting similar framework fragments.
    guest = fragment("acetylene.xyz", fragment_path=MOIETIES)
    matches = find_pattern(guest, repaired, disconnected_component=True)
    if not matches.nb_locations():
        raise SystemExit("Expected isolated acetylene guests were not found.")
    child = replace_pattern(matches, None, random_seed=args.seed)
    child.write_cif(out / "cleaned_SIFSIX.cif")
    write_report(
        out,
        parent_atoms=parent.n_atoms,
        repaired_atoms=repaired.n_atoms,
        child_atoms=child.n_atoms,
        disordered_locations=disorder.nb_locations(),
        removed_guests=matches.nb_locations(),
        validation=validate_structure(infer_bonds(child)).to_dict(),
    )


if __name__ == "__main__":
    main()
