#!/usr/bin/env python3
"""Compare symmetry after modification and construct a supercell.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, MOIETIES, example_parser, output_dir, write_report


def main():
    args = example_parser("symmetry_analysis", __doc__, seed=True).parse_args()
    out = output_dir(args)
    import numpy as np
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    from mofforge import Crystal, find_pattern, fragment, infer_bonds, replace_pattern

    parent = infer_bonds(Crystal.from_cif(CRYSTALS / "NiPyC_fragment_trouble.cif"))
    query = fragment("PyC.xyz", fragment_path=MOIETIES)
    replacement = fragment("PyC-CH3.xyz", fragment_path=MOIETIES)
    match = find_pattern(query, parent)
    if not match.nb_locations():
        raise SystemExit("Expected PyC linker was not found.")
    child = replace_pattern(match, replacement, nb_loc=1, random_seed=args.seed)
    # Wrapping changes coordinate representatives, preserving periodic positions.
    wrapped = child.wrap()
    np.testing.assert_allclose(
        child.lattice.get_all_distances(child.frac_coords, wrapped.frac_coords).diagonal(),
        0,
        atol=1e-8,
    )
    structure = wrapped.structure.copy()
    structure.make_supercell([2, 1, 1])
    supercell = infer_bonds(Crystal.from_structure(structure, name="NiPyC_supercell"))
    child.write_cif(out / "NiPyC_CH3.cif")
    supercell.write_cif(out / "NiPyC_CH3_supercell.cif")
    write_report(
        out,
        seed=args.seed,
        parent_atoms=parent.n_atoms,
        child_atoms=child.n_atoms,
        supercell_atoms=supercell.n_atoms,
        symprec_angstrom=0.1,
        parent_spacegroup=SpacegroupAnalyzer(
            parent.structure, symprec=0.1
        ).get_space_group_symbol(),
        child_spacegroup=SpacegroupAnalyzer(child.structure, symprec=0.1).get_space_group_symbol(),
        sample_fractional=wrapped.frac_coords[0],
        sample_cartesian_angstrom=wrapped.cart_coords[0],
    )


if __name__ == "__main__":
    main()
