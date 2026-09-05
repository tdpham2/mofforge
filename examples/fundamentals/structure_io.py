#!/usr/bin/env python3
"""Read CIF/XYZ data, infer bonds, and export structures.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, MOIETIES, example_parser, output_dir, write_report


def main():
    args = example_parser("structure_io", __doc__).parse_args()
    out = output_dir(args)
    from mofforge import Crystal, fragment, infer_bonds, read_xyz

    # CIF preserves the periodic cell; bond inference is a separate operation.
    parent = infer_bonds(Crystal.from_cif(CRYSTALS / "IRMOF-1.cif"), periodic=True)
    query = fragment("p-phenylene.xyz", fragment_path=MOIETIES)
    parent.write_cif(out / "framework.cif")
    query.write_xyz(out / "phenylene.xyz")
    species, coordinates = read_xyz(out / "phenylene.xyz")
    reloaded = Crystal.from_cif(out / "framework.cif")
    zinc = parent.species.index("Zn")
    write_report(
        out,
        atoms=parent.n_atoms,
        bonds=parent.n_bonds,
        formula=parent.structure.composition.formula,
        lattice_angstrom=parent.lattice.abc,
        zinc_coordination=parent.coordination_number(zinc),
        fragment_atoms=len(species),
        fragment_shape=coordinates.shape,
        reloaded_atoms=reloaded.n_atoms,
    )


if __name__ == "__main__":
    main()
