#!/usr/bin/env python3
"""Identify disconnected guests and compare solvent-retention controls.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, example_parser, output_dir, write_report


def main():
    args = example_parser("solvent_removal", __doc__, seed=True).parse_args()
    out = output_dir(args)
    from mofforge import Crystal, infer_bonds, place_adsorbate, remove_solvent

    parent = infer_bonds(Crystal.from_cif(CRYSTALS / "IRMOF-1.cif"))
    # Construct a known loaded input rather than relying on unidentified CIF guests.
    loaded = place_adsorbate(parent, "CO2", n_adsorbates=2, grid_spacing=1.0, random_seed=args.seed)
    loaded.crystal.write_cif(out / "loaded.cif")
    cleaned = remove_solvent(loaded.crystal)
    # min_atoms is a retention threshold: 3 keeps these three-atom molecules.
    retained = remove_solvent(loaded.crystal, min_atoms=3)
    explicit = remove_solvent(loaded.crystal, n_framework_components=1, keep_metal_containing=True)
    cleaned.crystal.write_cif(out / "desolvated.cif")
    write_report(
        out,
        parent_atoms=parent.n_atoms,
        placed=loaded.n_adsorbates,
        loaded_atoms=loaded.crystal.n_atoms,
        cleaned_atoms=cleaned.crystal.n_atoms,
        removed_atoms=cleaned.n_atoms_removed,
        removed_components=cleaned.n_components_removed,
        removed_formulas=[m.formula for m in cleaned.removed_molecules],
        threshold_retained_atoms=retained.crystal.n_atoms,
        explicit_framework_atoms=explicit.crystal.n_atoms,
    )


if __name__ == "__main__":
    main()
